import sys
import os
import pytest
import math

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from main import (
    read_sequence,
    bit_frequency_analysis,
    identical_consecutive_bits,
    longest_run_blocks_test,
    run_tests_on_files,
    print_results,
)

# read_sequence tests
def test_read_sequence_ok(tmp_path):
    path = tmp_path / "seq.txt"
    path.write_text("101010")
    assert read_sequence(str(path)) == "101010"

def test_read_sequence_file_not_found(monkeypatch):
    printed = []
    monkeypatch.setattr("builtins.print", lambda x: printed.append(x))
    result = read_sequence("no_such_file.txt")
    assert result == ""
    assert any("не найден" in s for s in printed)

def test_read_sequence_ioerror(monkeypatch):
    def fake_open(*args, **kwargs):
        raise IOError("Ошибка ввода-вывода")
    monkeypatch.setattr("builtins.open", fake_open)
    printed = []
    monkeypatch.setattr("builtins.print", lambda x: printed.append(x))
    result = read_sequence("file.txt")
    assert result == ""
    assert any("Ошибка при чтении файла" in s for s in printed)

# bit_frequency_analysis tests
@pytest.mark.parametrize(
    "seq, expected_sum",
    [("1",1),("0",-1),("10",0),("1100",0),("111000",0)]
)
def test_bit_frequency_sum(seq, expected_sum):
    length = len(seq)
    s_n = expected_sum / math.sqrt(length)
    expected = math.erfc(s_n / math.sqrt(2))
    result = bit_frequency_analysis(seq)
    assert pytest.approx(result, rel=1e-6) == expected

def test_bit_frequency_empty_string():
    with pytest.raises(ZeroDivisionError):
        bit_frequency_analysis("")

# identical_consecutive_bits tests
def test_identical_consecutive_bits_valid_case():
    seq = "10101010"
    result = identical_consecutive_bits(seq)
    assert result > 0

def test_identical_consecutive_bits_invalid_case():
    seq = "11111110000000"
    result = identical_consecutive_bits(seq)
    assert result < 0.01

# longest_run_blocks_test tests
@pytest.mark.parametrize(
    "seq",
    ["00000000","10000000","11000000","11100000","11110000"]
)
def test_longest_run_blocks_test_categories(seq):
    sequence = seq*50
    p = longest_run_blocks_test(sequence)
    assert 0 <= p <= 1

# run_tests_on_files
def test_run_tests_on_files_multi(monkeypatch):
    def fake_read(filename):
        if "a" in filename:
            return "1010"*10
        else:
            return "11110000"*10
    monkeypatch.setattr("main.read_sequence", fake_read)
    files = ["a.txt","b.txt"]
    results = run_tests_on_files(files)
    for name in ["a","b"]:
        assert name in results["frequency"]
        assert name in results["runs"]
        assert name in results["longest_run"]

# print_results
def test_print_results_output(capsys):
    results = {
        "frequency":{"test":0.123456},
        "runs":{"test":0.654321},
        "longest_run":{"test":0.02},
    }
    print_results(results)
    captured = capsys.readouterr().out
    assert "Частотный побитовый тест" in captured
    assert "0.123456" in captured
    assert "0.654321" in captured
    assert "0.020000" in captured
    assert "Случайна" in captured

# integration tests with real sequences
def test_run_tests_on_real_files(tmp_path):
    seq_cpp = "00001100110100001111111110111001100000010100101111000010001101000000111111011101110101111000010100110011010100100011100111000111"
    seq_java = "11001110011100110000000111011110110101101100000101111001100001100110111101111010101111010011110111000011001010010001001000011000"

    file_cpp = tmp_path / "sequence_cpp.txt"
    file_java = tmp_path / "sequence_java.txt"
    file_cpp.write_text(seq_cpp)
    file_java.write_text(seq_java)

    files = [str(file_cpp), str(file_java)]
    results = run_tests_on_files(files)

    for name in ["sequence_cpp", "sequence_java"]:
        assert name in results["frequency"]
        assert name in results["runs"]
        assert name in results["longest_run"]

    for category in ["frequency","runs","longest_run"]:
        for value in results[category].values():
            assert 0 <= value <= 1
