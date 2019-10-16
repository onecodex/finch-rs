use std::fs::File;
use std::io::BufReader;

use numpy::{PyArray, PyArray1, PyArray2};
use pyo3::class::*;
use pyo3::prelude::*;
use pyo3::types::{PyBytes, PyTuple, PyType};
use pyo3::{py_exception, wrap_function};

use crate::distance::{common_counts, distance, distance_scaled, minmer_matrix};
use crate::filtering::FilterParams;
use crate::serialization::{
    read_finch_file, read_mash_file, write_finch_file, FinchSketch, MultiSketch as MSType,
    FINCH_BIN_EXT, MASH_EXT,
};
use crate::sketch_files as rs_sketch_files;
use crate::sketch_schemes::{KmerCount, SketchParams};
use crate::Result as FinchResult;

py_exception!(finch, FinchError, pyo3::exceptions::Exception);

#[pyclass]
/// A Multisketch is a collection of Sketchs with information about their generation parameters
/// (to make sure they're consistant for distance calculation).
pub struct Multisketch {
    pub ms: MSType,
}

#[pymethods]
impl Multisketch {
    #[classmethod]
    /// open_file(filename)
    ///
    /// Takes a file path to either a `.sk` or a `.mash` file and returns the Multisketch
    /// represented by that file.
    fn open_file(_cls: &PyType, filename: &str) -> PyResult<Multisketch> {
        // this is the same as in main.rs and we should probably refactor these together
        let file = File::open(filename)?;
        let mut buf_reader = BufReader::new(file);
        let ms: MSType = if filename.ends_with(MASH_EXT) {
            read_mash_file(&mut buf_reader)
                .map_err(|e| PyErr::new::<FinchError, _>(format!("{}", e)))?
        } else if filename.ends_with(FINCH_BIN_EXT) {
            read_finch_file(&mut buf_reader)
                .map_err(|e| PyErr::new::<FinchError, _>(format!("{}", e)))?
        } else {
            // TODO: check for a finch extension?
            serde_json::from_reader(buf_reader)
                .map_err(|e| PyErr::new::<FinchError, _>(format!("{}", e)))?
        };
        Ok(Multisketch { ms: ms })
    }

    // TODO: save method?

    #[getter]
    fn get_kmer(&self) -> PyResult<u8> {
        Ok(self.ms.kmer)
    }

    #[getter]
    fn get_alphabet(&self) -> PyResult<String> {
        Ok(self.ms.alphabet.clone())
    }

    #[getter]
    fn get_preserve_case(&self) -> PyResult<bool> {
        Ok(self.ms.preserve_case)
    }

    #[getter]
    fn get_canonical(&self) -> PyResult<bool> {
        Ok(self.ms.canonical)
    }

    #[getter]
    fn get_sketch_size(&self) -> PyResult<u32> {
        Ok(self.ms.sketch_size)
    }

    #[getter]
    fn get_hash_type(&self) -> PyResult<String> {
        Ok(self.ms.hash_type.clone())
    }

    #[getter]
    fn get_hash_bits(&self) -> PyResult<u16> {
        Ok(self.ms.hash_bits)
    }

    #[getter]
    fn get_hash_seed(&self) -> PyResult<u64> {
        Ok(self.ms.hash_seed)
    }

    #[getter]
    fn get_sketches(&self) -> PyResult<Vec<Sketch>> {
        // TODO: we should be doing this without a clone probably?
        // (i.e. returning references to things in the multisketch)
        let sketches: Vec<FinchSketch> = (&self.ms).into();
        let fsketches: FinchResult<Vec<Sketch>> =
            sketches.into_iter().map(|s| Ok(Sketch { s })).collect();

        fsketches.map_err(|e| PyErr::new::<FinchError, _>(format!("{}", e)))
    }
}

#[pyproto]
impl PyObjectProtocol for Multisketch {
    fn __repr__(&self) -> PyResult<String> {
        let n_sketches = self.ms.sketches.len();
        let sketch_plural = if n_sketches == 1 {
            "sketch"
        } else {
            "sketches"
        };
        Ok(format!("<Multisketch ({} {})>", n_sketches, sketch_plural))
    }
}

#[pyproto]
impl PyMappingProtocol for Multisketch {
    fn __len__(&self) -> PyResult<usize> {
        Ok(self.ms.sketches.len())
    }
}

/// A Sketch is a collection of deterministically-selected hashes from a single
/// sequencing file.
#[pyclass]
pub struct Sketch {
    pub s: FinchSketch,
}

#[pymethods]
impl Sketch {
    #[new]
    fn __new__(obj: &PyRawObject, name: &str) -> PyResult<()> {
        // TODO: take a hashes parameter: Vec<(usize, &[u8], u16, u16)>,
        // and a sketch_params?
        let sketch_params = SketchParams::Mash {
            kmers_to_sketch: 1000,
            final_size: 1000,
            no_strict: true,
            kmer_length: 21,
            hash_seed: 0,
        };
        let sketch = FinchSketch {
            name: name.to_string(),
            seq_length: 0,
            num_valid_kmers: 0,
            comment: String::new(),
            hashes: Vec::new(),
            sketch_params,
            filter_params: FilterParams::default(),
        };
        obj.init(|_| Sketch { s: sketch })
    }

    #[getter]
    fn get_name(&self) -> PyResult<String> {
        Ok(self.s.name.clone())
    }

    #[setter]
    fn set_name(&mut self, value: &str) -> PyResult<()> {
        self.s.name = value.to_string();
        Ok(())
    }

    #[getter]
    fn get_seq_length(&self) -> PyResult<u64> {
        Ok(self.s.seq_length)
    }

    #[getter]
    fn get_num_valid_kmers(&self) -> PyResult<u64> {
        Ok(self.s.num_valid_kmers)
    }

    #[getter]
    fn get_comment(&self) -> PyResult<String> {
        Ok(self.s.comment.clone())
    }

    #[getter]
    fn get_hashes(&self) -> PyResult<Vec<(u64, Py<PyBytes>, u64, u64)>> {
        let gil = Python::acquire_gil();
        let py = gil.python();
        Ok(self
            .s
            .hashes
            .clone()
            .into_iter()
            .map(|i| (i.hash, PyBytes::new(py, &i.kmer), i.count, i.extra_count))
            .collect())
    }

    // TODO: there are a lot of issues to fix in here; we should also try to destructure the
    // list depending on the format of the tuples; i.e. allow (usize, &[u8], u16, u16),
    // (usize, &[u8], u16), (usize, u16), (usize, &[u8]), usize, etc

    // #[setter]
    // fn set_hashes(&self, value: &PyObject) -> PyResult<()> {

    //     let value: &[(usize, PyBytes, u16, u16)] = PyObjectRef::extract(value)?;
    //     let kmers: Vec<KmerCount> = value.iter().map(|(hash, kmer, count, extra_count)| {
    //         KmerCount {
    //             hash: *hash,
    //             kmer: kmer.as_bytes().to_vec(),
    //             count: *count,
    //             extra_count: *extra_count,
    //         }
    //     }).collect();
    //     self.s.set_kmers(&kmers);
    //     Ok(())
    // }

    // TODO: filtering method

    // TODO: clip to n kmers/hashes method

    // /// merge(sketch, size, scale)
    // ///
    // /// Merge the second sketch into this one. If size is specified, use
    // /// that as the new sketch's size. If scale is specified, merge the
    // /// sketches together as if they are scaled sketches (for scaled sketches
    // /// that have "high" hashes because they're under `size`, this will
    // /// potentially remove those hashes if the new sketch is large enough).
    pub fn merge(&mut self, sketch: &Sketch, size: Option<usize>) -> PyResult<()> {
        // update my parameters from the remote's
        self.s.seq_length += sketch.s.seq_length;
        self.s.num_valid_kmers += sketch.s.num_valid_kmers;

        // TODO: do something with filters?
        // TODO: we should also check the sketch_params are compatible?

        // now merge the hashes together; someday it would be nice to use something idiomatic like:
        // https://users.rust-lang.org/t/solved-merge-multiple-sorted-vectors-using-iterators/6543
        let sketch1 = &self.s.hashes;
        let sketch2 = &sketch.s.hashes;

        let mut new_hashes = Vec::with_capacity(sketch1.len() + sketch2.len());
        let (mut i, mut j) = (0, 0);
        while (i < sketch1.len()) && (j < sketch2.len()) {
            if sketch1[i].hash < sketch2[j].hash {
                new_hashes.push(sketch1[i].clone());
                i += 1;
            } else if sketch2[j].hash < sketch1[i].hash {
                new_hashes.push(sketch2[j].clone());
                j += 1;
            } else {
                new_hashes.push(KmerCount {
                    hash: sketch1[i].hash,
                    kmer: sketch1[i].kmer.clone(),
                    count: sketch1[i].count + sketch2[j].count,
                    extra_count: sketch1[i].extra_count + sketch2[j].extra_count,
                });
                i += 1;
                j += 1;
            }
        }

        // now clip to the appropriate size
        let scale = self.s.sketch_params.hash_info().3;
        match (size, scale) {
            (Some(s), Some(sc)) => {
                let max_hash = u64::max_value() / (1. / sc) as u64;
                // truncate to hashes <= max/sc (or) s whichever is higher
                new_hashes = new_hashes
                    .iter()
                    .enumerate()
                    .take_while(|(ix, h)| (h.hash <= max_hash) || (*ix < s))
                    .map(|(_, h)| h.clone())
                    .collect();
            }
            (None, Some(sc)) => {
                let max_hash = u64::max_value() / (1. / sc) as u64;
                // truncate to hashes <= max/sc
                new_hashes = new_hashes
                    .iter()
                    .take_while(|h| h.hash <= max_hash)
                    .map(|h| h.clone())
                    .collect();
            }
            (Some(s), None) => {
                // truncate to size
                new_hashes.truncate(s);
            }
            (None, None) => {
                // no filtering
            }
        }
        self.s.hashes = new_hashes;
        Ok(())
    }

    /// compare(sketch, mash_mode=False)
    ///
    /// Calculates the containment within and jaccard similarity to another sketch.
    #[args(mash_mode = true)]
    pub fn compare(&self, sketch: &Sketch, mash_mode: bool) -> PyResult<(f64, f64)> {
        let dist = distance(&self.s.hashes, &sketch.s.hashes, &"", &"", mash_mode)
            .map_err(|e| PyErr::new::<FinchError, _>(format!("{}", e)))?;

        Ok((dist.containment, dist.jaccard))
    }

    /// compare_scaled(sketch)
    ///
    /// Calculates the containment within and jaccard similarity to another scaled sketch.
    pub fn compare_scaled(&self, sketch: &Sketch) -> PyResult<(f64, f64)> {
        let dist = distance_scaled(&self.s.hashes, &sketch.s.hashes, &"", &"")
            .map_err(|e| PyErr::new::<FinchError, _>(format!("{}", e)))?;

        Ok((dist.containment, dist.jaccard))
    }

    /// compare_counts(sketch)
    ///
    ///
    pub fn compare_counts(&self, sketch: &Sketch) -> PyResult<(u64, u64, u64, u64, u64)> {
        Ok(common_counts(&self.s.hashes, &sketch.s.hashes))
    }

    // /// compare_matrix(*sketches)
    // ///
    // /// Generates a numpy matrix of hash/kmer counts aligned to a "primary"
    // /// reference. This matrix can then be used for downstream NNLS analysis.
    #[args(args = "*")]
    pub fn compare_matrix(&self, args: &PyTuple) -> PyResult<Py<PyArray2<u64>>> {
        let sketches: Vec<&Sketch> = args.extract()?;
        let sketch_kmers: Vec<&[KmerCount]> = sketches.iter().map(|s| &s.s.hashes[..]).collect();
        let result = minmer_matrix(&self.s.hashes, &sketch_kmers);

        let gil = Python::acquire_gil();
        let py = gil.python();
        Ok(PyArray::from_owned_array(py, result).to_owned())
    }

    #[getter]
    pub fn get_counts(&self) -> PyResult<Py<PyArray1<u64>>> {
        let result = self.s.hashes.iter().map(|k| u64::from(k.count));

        let gil = Python::acquire_gil();
        let py = gil.python();
        Ok(PyArray::from_exact_iter(py, result).to_owned())
    }
}

#[pyproto]
impl PyObjectProtocol for Sketch {
    fn __repr__(&self) -> PyResult<String> {
        Ok(format!("<Sketch \"{}\">", self.s.name.clone()))
    }
}

#[pyproto]
impl PyMappingProtocol for Sketch {
    fn __len__(&self) -> PyResult<usize> {
        Ok(self.s.len())
    }
}

// TODO: impl PyNumberProtocol addition or subtraction for Sketch to allow merging/
// set difference calculations?
// see https://github.com/PyO3/pyo3/blob/master/tests/test_arithmetics.rs for details

/// sketch_files(filenames, n_hashes, final_size, kmer_length, filter, seed)
/// ---
///
/// From the FASTA and FASTQ file paths, create a Multisketch.
// #[pyfunction(n_hashes=null, kmer_length=21, filter=true, seed=0)]  // TODO: this doesn't work?
#[pyfunction]
pub fn sketch_files(
    filenames: Vec<&str>,
    n_hashes: usize,
    final_size: Option<usize>,
    kmer_length: u8,
    filter: bool,
    seed: u64,
) -> PyResult<Multisketch> {
    // TODO: allow more filter customization?

    // TODO: allow passing in a single file without the list
    let sketch_params = SketchParams::Mash {
        kmers_to_sketch: n_hashes,
        final_size: final_size.unwrap_or(n_hashes),
        no_strict: false,
        kmer_length,
        hash_seed: seed,
    };
    let filters = FilterParams {
        filter_on: Some(filter),
        abun_filter: (None, None),
        err_filter: 1.,
        strand_filter: 0.1,
    };

    Ok(Multisketch {
        ms: rs_sketch_files(&filenames, &sketch_params, &filters)
            .map_err(|e| PyErr::new::<FinchError, _>(format!("{}", e)))?,
    })
}

#[pyfunction]
pub fn write_sketch_file(filename: &str, sketches: Vec<&Sketch>) -> PyResult<()> {
    let fsketches: Vec<FinchSketch> = sketches.iter().map(|s| s.s.clone()).collect();

    let mut out = File::create(&filename)
        .map_err(|_| PyErr::new::<FinchError, _>(format!("Could not create {}", filename)))?;
    write_finch_file(&mut out, &fsketches)
        .map_err(|e| PyErr::new::<FinchError, _>(format!("{}", e)))?;
    Ok(())
}

#[pymodinit]
/// Finch is a MinHash sketch processing library.
fn finch(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<Multisketch>()?;
    m.add_class::<Sketch>()?;
    m.add_function(wrap_function!(sketch_files))?;
    m.add_function(wrap_function!(write_sketch_file))?;

    Ok(())
}
