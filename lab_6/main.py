import math
import argparse
from scipy.special import gammaincc
import os

def read_sequence(filename: str) -> str:
    try:
        with open(filename, 'r') as file:
            sequence = file.read().strip()
        return sequence
    except FileNotFoundError:
        print(f"Ошибка: файл '{filename}' не найден.")
        return ""
    except IOError:
        print(f"Ошибка при чтении файла '{filename}'.")
        return ""

def bit_frequency_analysis(sequence: str) -> float:
    sequence_sum = 0
    length = len(sequence)
    for i in sequence:
        match i:
            case '1':
                sequence_sum += 1
            case '0':
                sequence_sum -= 1
    s_n = (1 / math.sqrt(length)) * sequence_sum
    p_value = math.erfc(s_n / math.sqrt(2))
    return p_value

def identical_consecutive_bits(sequence: str) -> float:
    length = len(sequence)
    share_of_ones = sequence.count('1') / length
    if abs(share_of_ones - 0.5) < (2 / math.sqrt(length)):
        v_n = 0
        for i in range(length - 1):
            if sequence[i] != sequence[i + 1]:
                v_n += 1
        numerator = abs(v_n - 2 * length * share_of_ones * (1 - share_of_ones))
        denominator = 2 * math.sqrt(2 * length) * share_of_ones * (1 - share_of_ones)
        p_value = math.erfc(numerator / denominator)
        return p_value
    else:
        return 0

def longest_run_blocks_test(sequence: str, block_size: int = 8) -> float:
    pi = [0.2148, 0.3672, 0.2305, 0.1875]
    blocks = [sequence[i:i + block_size] for i in range(0, len(sequence), block_size)]
    v = [0, 0, 0, 0]
    for block in blocks:
        max_run = 0
        current_run = 0
        for bit in block:
            match bit:
                case '1':
                    current_run += 1
                    max_run = max(max_run, current_run)
                case '0':
                    current_run = 0
        match max_run:
            case r if r <= 1:
                v[0] += 1
            case 2:
                v[1] += 1
            case 3:
                v[2] += 1
            case r if r >= 4:
                v[3] += 1
    chi_square = sum((v[i] - len(blocks) * pi[i]) ** 2 / (len(blocks) * pi[i]) for i in range(4))
    p_value = gammaincc(1.5, chi_square / 2)
    return p_value

def run_tests_on_files(files):
    results = {'frequency': {}, 'runs': {}, 'longest_run': {}}
    for filename in files:
        sequence = read_sequence(filename)
        base_name = os.path.splitext(os.path.basename(filename))[0]
        results['frequency'][base_name] = bit_frequency_analysis(sequence)
        results['runs'][base_name] = identical_consecutive_bits(sequence)
        results['longest_run'][base_name] = longest_run_blocks_test(sequence)
    return results

def print_results(results):
    print("Частотный побитовый тест:")
    for name, p_value in results['frequency'].items():
        print(f"{name}: {p_value:.6f}")

    print("\nТест на одинаковые подряд идущие биты:")
    for name, p_value in results['runs'].items():
        print(f"{name}: {p_value:.6f}")

    print("\nТест на самую длинную последовательность единиц в блоке (size=8):")
    for name, p_value in results['longest_run'].items():
        print(f"{name}: {p_value:.6f}")
        print(f"Заключение: {'Случайна' if p_value >= 0.01 else 'Неслучайна'}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Тестирование последовательностей по стандарту NIST.')
    parser.add_argument('files', nargs='+', help='Пути к файлам с последовательностями')
    args = parser.parse_args()
    results = run_tests_on_files(args.files)
    print_results(results)
