use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;

use numpy::{PyArray, PyArray1, PyArray2};
use pyo3::class::*;
use pyo3::prelude::*;
use pyo3::types::{PyBytes, PyTuple, PyType};
use pyo3::{py_exception, wrap_function};

use crate::distance::{common_counts, distance, distance_scaled, minmer_matrix};
use crate::filtering::FilterParams;
use crate::hash_schemes::KmerCount;
use crate::mash_files;
use crate::serialization::{read_mash_file, MultiSketch as MSType, Sketch as SType, MASH_EXT};

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
        } else {
            // TODO: check for a finch extension?
            serde_json::from_reader(buf_reader)
                .map_err(|e| PyErr::new::<FinchError, _>(format!("{}", e)))?
        };
        Ok(Multisketch { ms: ms })
    }

    // TODO: save method

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
        Ok(self.ms.preserveCase)
    }

    #[getter]
    fn get_canonical(&self) -> PyResult<bool> {
        Ok(self.ms.canonical)
    }

    #[getter]
    fn get_sketch_size(&self) -> PyResult<u32> {
        Ok(self.ms.sketchSize)
    }

    #[getter]
    fn get_hash_type(&self) -> PyResult<String> {
        Ok(self.ms.hashType.clone())
    }

    #[getter]
    fn get_hash_bits(&self) -> PyResult<u16> {
        Ok(self.ms.hashBits)
    }

    #[getter]
    fn get_hash_seed(&self) -> PyResult<u64> {
        Ok(self.ms.hashSeed)
    }

    #[getter]
    fn get_sketches(&self) -> PyResult<Vec<Sketch>> {
        // TODO: we should be doing this without a clone probably?
        Ok(self
            .ms
            .sketches
            .iter()
            .map(|s| Sketch { s: s.clone() })
            .collect())
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
    pub s: SType,
}

#[pymethods]
impl Sketch {
    #[new]
    #[args(length = 0, n_kmers = 0)]
    fn __new__(obj: &PyRawObject, name: &str, length: u64, n_kmers: u64) -> PyResult<()> {
        // TODO: take a hashes parameter: Vec<(usize, &[u8], u16, u16)>,
        let sketch = SType::new(name, length, n_kmers, Vec::new(), &HashMap::new());
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
    fn get_seq_length(&self) -> PyResult<Option<u64>> {
        Ok(self.s.seqLength)
    }

    #[getter]
    fn get_num_valid_kmers(&self) -> PyResult<Option<u64>> {
        Ok(self.s.numValidKmers)
    }

    #[getter]
    fn get_comment(&self) -> PyResult<Option<String>> {
        Ok(self.s.comment.clone())
    }

    #[getter]
    fn get_hashes(&self) -> PyResult<Vec<(u64, Py<PyBytes>, u16, u16)>> {
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

    // TODO: merge sketches method

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

/// sketch_files(filenames, n_hashes, final_size, kmer_length, filter, seed, scaled)
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
    scaled: Option<f64>,
) -> PyResult<Multisketch> {
    // TODO: allow more filter customization?

    // TODO: allow passing in a single file without the list
    let mut filters = FilterParams {
        filter_on: Some(filter),
        abun_filter: (None, None),
        err_filter: 1.,
        strand_filter: 0.1,
    };
    let fsize = match final_size {
        Some(x) => x,
        None => n_hashes,
    };

    Ok(Multisketch {
        ms: mash_files(
            &filenames,
            n_hashes,
            fsize,
            kmer_length,
            &mut filters,
            false,
            seed,
            scaled,
        )
        .map_err(|e| PyErr::new::<FinchError, _>(format!("{}", e)))?,
    })
}

#[pymodinit]
/// Finch is a MinHash sketch processing library.
fn finch(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<Multisketch>()?;
    m.add_class::<Sketch>()?;
    m.add_function(wrap_function!(sketch_files))?;

    Ok(())
}
