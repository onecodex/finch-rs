from pathlib import Path

import numpy as np
import pytest
from finch import FinchError, Multisketch, Sketch, sketch_file

QUERY_FILE = Path(__file__).resolve().parent / "cli/tests/data/query.fa"
REFS_FILE = Path(__file__).resolve().parent / "cli/tests/data/refs.fa"


@pytest.fixture
def small_sketch():
    # query.fa is small and so no_strict needs to be enabled otherwise the sequences will get
    # rejected as "too few k-mers"
    return sketch_file(QUERY_FILE.as_posix(), no_strict=True)


@pytest.fixture
def small_named_sketches():
    sketches = []

    for name in ("a", "b", "c"):
        s = sketch_file(QUERY_FILE.as_posix(), no_strict=True)
        s.name = name
        sketches.append(s)

    return sketches


@pytest.fixture
def sketch():
    return sketch_file(REFS_FILE.as_posix(), filter=False)


#
## sketch_file tests
#


def test_sketch_file_no_strict():
    s = sketch_file(QUERY_FILE.as_posix(), no_strict=True)

    assert isinstance(s, Sketch)
    assert len(s) > 0
    assert s.seq_length > 0
    assert s.num_valid_kmers > 0
    assert s.name == QUERY_FILE.as_posix()
    assert s.sketch_params["kmer_length"] == 21


@pytest.mark.parametrize(
    "input_path,kmer_length,n_hashes",
    [(QUERY_FILE, 11, 5), (REFS_FILE, 11, 5), (REFS_FILE, 21, 10)],
)
def test_sketch_file_n_hashes(input_path, kmer_length, n_hashes):
    s = sketch_file(
        input_path.as_posix(), n_hashes=n_hashes, kmer_length=kmer_length, filter=False
    )

    assert len(s) <= n_hashes
    assert len(s.hashes) <= n_hashes


@pytest.mark.parametrize(
    "input_path,kmer_length",
    [(REFS_FILE, 11), (REFS_FILE, 21), (REFS_FILE, 31)],
)
def test_sketch_file_kmer_length(input_path, kmer_length):
    s = sketch_file(input_path.as_posix(), kmer_length=kmer_length, filter=False)

    assert s.sketch_params["kmer_length"] == kmer_length


@pytest.mark.parametrize(
    "input_path,kmer_length",
    [(QUERY_FILE, 11), (QUERY_FILE, 7), (QUERY_FILE, 17)],
)
def test_sketch_file_kmer_length_no_strict(input_path, kmer_length):
    s = sketch_file(
        input_path.as_posix(), kmer_length=kmer_length, filter=False, no_strict=True
    )

    assert s.sketch_params["kmer_length"] == kmer_length


def test_sketch_file_missing_file_raises():
    with pytest.raises(FinchError):
        sketch_file("/does/not/exist.fa")


def test_sketch_file_too_few_kmers_raises():
    # With the default (strict) parameters the tiny test file is rejected.
    with pytest.raises(FinchError):
        sketch_file(QUERY_FILE.as_posix())


#
## Sketch struct tests
#


def test_new_sketch_is_empty():
    s = Sketch("empty")
    assert s.name == "empty"
    assert len(s) == 0
    assert s.seq_length == 0
    assert s.num_valid_kmers == 0
    assert s.hashes == []


def test_sketch_repr():
    assert repr(Sketch("foo")) == '<Sketch "foo">'


def test_sketch_name_setter():
    s = Sketch("foo")
    s.name = "bar"

    assert s.name == "bar"


def test_sketch_comment_setter():
    s = Sketch("foo")
    s.comment = "a helpful comment"

    assert s.comment == "a helpful comment"


def test_sketch_params(small_sketch):
    params = small_sketch.sketch_params

    assert params["sketch_type"] == "mash"
    assert params["kmer_length"] == 21
    assert "hash_seed" in params


def test_sketch_copy(small_sketch):
    clone = small_sketch.copy()
    clone.name = "clone"

    assert len(clone) == len(small_sketch)
    assert small_sketch.name != "clone"


@pytest.mark.parametrize("sketch_fix", ["small_sketch", "sketch"])
def test_sketch_compare(sketch_fix, request):
    sketch = request.getfixturevalue(sketch_fix)
    containment, jaccard = sketch.compare(sketch)

    assert containment == 1.0
    assert jaccard == 1.0


def test_sketch_compare_bounds(small_named_sketches):
    a, b, c = small_named_sketches
    containment, jaccard = a.compare(b)

    assert 0.0 <= containment <= 1.0
    assert 0.0 <= jaccard <= 1.0

    containment, jaccard = a.compare(c)

    assert 0.0 <= containment <= 1.0
    assert 0.0 <= jaccard <= 1.0


def test_sketch_compare_counts(small_sketch):
    result = small_sketch.compare_counts(small_sketch)

    assert len(result) == 8

    common = result[0]

    assert common == len(small_sketch)


@pytest.mark.parametrize("sketch_fix", ["small_sketch", "sketch"])
def test_compare_matrix_shape(sketch_fix, request):
    sketch = request.getfixturevalue(sketch_fix)
    other = sketch.copy()
    matrix = sketch.compare_matrix(other, other)

    assert isinstance(matrix, np.ndarray)
    assert matrix.dtype == np.int32

    # one row per query sketch, one column per reference hash
    assert matrix.shape == (2, len(sketch))


def test_merge_identical_keeps_size(small_sketch):
    other = small_sketch.copy()
    n_before = len(small_sketch)
    small_sketch.merge(other, None)
    # merging identical sketches dedupes the hashes, so the size is unchanged
    assert len(small_sketch) == n_before


def test_merge_sums_shared_counts(small_sketch):
    before = small_sketch.counts.copy()
    small_sketch.merge(small_sketch.copy(), None)
    assert (small_sketch.counts == before * 2).all()


#
## Multisketch struct tests
#


def test_multisketch_from_sketches(small_named_sketches):
    ms = Multisketch.from_sketches(small_named_sketches)

    assert len(ms) == 3
    assert repr(ms) == "<Multisketch (3 sketches)>"


def test_multisketch_getitem_by_index(small_named_sketches):
    ms = Multisketch.from_sketches(small_named_sketches)

    assert ms[0].name == "a"
    assert ms[2].name == "c"


def test_multisketch_getitem_by_name(small_named_sketches):
    ms = Multisketch.from_sketches(small_named_sketches)

    assert ms["b"].name == "b"


def test_multisketch_getitem_missing_name_raises(small_named_sketches):
    ms = Multisketch.from_sketches(small_named_sketches)

    with pytest.raises(KeyError):
        ms["does-not-exist"]


def test_multisketch_getitem_out_of_range_raises(small_named_sketches):
    ms = Multisketch.from_sketches(small_named_sketches)

    with pytest.raises(IndexError):
        ms[42]


def test_multisketch_contains(small_named_sketches):
    ms = Multisketch.from_sketches(small_named_sketches)

    assert "a" in ms
    assert "does-not-exist" not in ms


def test_multisketch_add(small_named_sketches):
    ms = Multisketch.from_sketches(small_named_sketches[:1])
    ms.add(small_named_sketches[1])

    assert len(ms) == 2
    assert ms[0].name == "a"
    assert ms[1].name == "b"


def test_multisketch_delitem(small_named_sketches):
    ms = Multisketch.from_sketches(small_named_sketches)

    del ms[0]

    assert len(ms) == 2
    assert [s.name for s in ms] == ["b", "c"]


def test_multisketch_best_match(small_named_sketches):
    ms = Multisketch.from_sketches(small_named_sketches)
    idx, match = ms.best_match(small_named_sketches[0])

    assert isinstance(idx, int)
    assert isinstance(match, Sketch)


def test_multisketch_filter_to_names(small_named_sketches):
    ms = Multisketch.from_sketches(small_named_sketches)
    ms.filter_to_names(["a", "c"])

    assert sorted(s.name for s in ms) == ["a", "c"]


def test_multisketch_filter_to_matches(small_named_sketches):
    ms = Multisketch.from_sketches(small_named_sketches)
    # identical sketches all match the query perfectly
    ms.filter_to_matches(small_named_sketches[0], 1.0)

    assert len(ms) == 3


def test_multisketch_save(small_named_sketches, tmp_path):
    ms = Multisketch.from_sketches(small_named_sketches)
    path = tmp_path / "sketches.bsk"
    ms.save(str(path))

    assert path.exists()

    reopened = Multisketch.open(str(path))

    assert len(reopened) == len(ms)
    assert sorted(s.name for s in reopened) == ["a", "b", "c"]


def test_multisketch_open_missing_raises():
    with pytest.raises(FinchError):
        Multisketch.open("/does/not/exist.bsk")
