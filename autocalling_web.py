from collections import defaultdict
import pandas as pd
import numpy as np
from Bio import SeqIO
from scipy.signal import find_peaks
from scipy.ndimage import uniform_filter1d
import glob, os
import matplotlib.pyplot as plt
import itertools

# ====================================================데이터 정의=====================================================
WINDOW_SIZE = 31
COLORS = ['blue', 'green', 'black', 'red']
CHANNEL_KEYS = ['DATA9', 'DATA10', 'DATA11', 'DATA12']

all_best_params = {
    'P1_blue': {'k_blue': 1.392833480861528, 'c_blue': -219},
    'P1_green': {'k_green': 1.6166427202524143, 'c_green': 619},
    'P1_black': {'k_black': 0.21911473524003966, 'c_black': 626},
    'P1_red': {'k_red': 0.44764811722857727, 'c_red': 333},
    'P2_blue': {'k_blue': 1.7286536958531877, 'c_blue': -150},
    'P2_green': {'k_green': 0.8745067389785239, 'c_green': 756},
    'P2_black': {'k_black': 1.11806008203225, 'c_black': 225},
    'P2_red': {'k_red': 1.3672599788222723, 'c_red': 620}
    }

alleles_ref = {'*1': ['G', 'C', 'T', 'G', 'G', 'G', 'G', 'C', 'G', 'A', 'G', 'A', 'C', 'A', 'T', 'G', 'G'],
               '*1XN': ['G', 'C', 'T', 'G', 'G', 'G', 'G', 'C', 'G', 'A', 'G', 'A', 'C', 'A', 'T', 'G', 'A'],
               '*2': ['G', 'C', 'T', 'G', 'G', 'G', 'G', 'T', 'G', 'A', 'G', 'A', 'C', 'A', 'T', 'G', 'G'],
               '*2XN': ['G', 'C', 'T', 'G', 'G', 'G', 'G', 'T', 'G', 'A', 'G', 'A', 'C', 'A', 'T', 'G', 'A'],
               '*3': ['G', 'C', 'T', 'G', 'G', 'G', 'G', 'C', 'G', 'A', 'G', 'G', 'C', 'A', 'T', 'G', 'G'],
               '*4': ['G', 'T', 'T', 'G', 'G', 'G', 'G', 'C', 'G', 'A', 'A', 'A', 'C', 'A', 'T', 'G', 'G'],
               '*4X2': ['G', 'T', 'T', 'G', 'G', 'G', 'G', 'C', 'G', 'G', 'A', 'A', 'C', 'A', 'T', 'G', 'A'],
               '*5': ['-', '-', '-', '-', '-', '-', '-', '-', '-', 'G', '-', '-', '-', '-', '-', '-', '-'],
               '*10D': ['G', 'T', 'T', 'G', 'G', 'G', 'G', 'C', 'G', '-', 'G', 'A', 'C', 'A', 'T', 'G', 'G'],
               '*6': ['G', 'C', 'T', 'G', 'G', 'G', 'G', 'C', 'G', 'A', 'G', 'A', 'C', 'A', 'G', 'G', 'G'],
               '*9': ['G', 'C', 'T', 'G', 'G', 'G', 'G', 'C', 'G', 'A', 'G', 'A', 'C', 'G', 'T', 'G', 'G'],
               '*10B': ['G', 'T', 'T', 'G', 'G', 'G', 'G', 'C', 'G', 'A', 'G', 'A', 'C', 'A', 'T', 'G', 'G'],
               '*10BX2': ['G', 'T', 'T', 'G', 'G', 'G', 'G', 'C', 'G', 'A', 'G', 'A', 'C', 'A', 'T', 'G', 'A'],
               '*14B': ['A', 'C', 'T', 'G', 'G', 'G', 'G', 'T', 'G', 'A', 'G', 'A', 'C', 'A', 'T', 'G', 'G'],
               '*17': ['G', 'C', 'T', 'G', 'G', 'G', 'G', 'T', 'G', 'A', 'G', 'A', 'T', 'A', 'T', 'G', 'G'],
               '*17XN': ['G', 'C', 'T', 'G', 'G', 'G', 'G', 'T', 'G', 'A', 'G', 'A', 'T', 'A', 'T', 'G', 'A'],
               '*18': ['G', 'C', 'T', 'G', 'G', 'G', 'T', 'C', 'G', 'A', 'G', 'A', 'C', 'A', 'T', 'G', 'G'],
               '*21B': ['G', 'C', 'T', 'C', 'G', 'G', 'G', 'T', 'G', 'A', 'G', 'A', 'C', 'A', 'T', 'G', 'G'],
               '*29': ['G', 'C', 'T', 'G', 'G', 'G', 'G', 'T', 'G', 'A', 'G', 'A', 'C', 'A', 'T', 'A', 'G'],
               '*41': ['G', 'C', 'T', 'G', 'A', 'G', 'G', 'T', 'G', 'A', 'G', 'A', 'C', 'A', 'T', 'G', 'G'],
               '*49': ['G', 'T', 'A', 'G', 'G', 'G', 'G', 'C', 'G', 'A', 'G', 'A', 'C', 'A', 'T', 'G', 'G'],
               '*52': ['G', 'T', 'T', 'G', 'G', 'A', 'G', 'C', 'G', 'A', 'G', 'A', 'C', 'A', 'T', 'G', 'G'],
               '*60': ['G', 'C', 'T', 'G', 'G', 'G', 'G', 'C', 'A', 'A', 'G', 'A', 'C', 'A', 'T', 'G', 'G']
               }

mutation_map = pd.DataFrame([
    ["1758G>A", "blue", "green", "b"],
    ["100C>T", "black", "red", "b"],
    ["1611T>A", "green", "red", "b"],
    ["2573_2574insC", "blue", "black", "b"],
    ["2988G>A", "blue", "green", "b"],
    ["3877G>A", "black", "red", "b"],
    ["4125_4133\ndupGTGCCCACT", "blue", "red", "b"],
    ["2850C>T", "blue", "green", "b"],
    ["1887insTA", "black", "red", "b"],
    ["Deletion", "red", "black", "f"],
    ["1846G>A", "blue", "green", "b"],
    ["2549delA", "green", "blue", "f"],
    ["1023C>T", "black", "red", "b"],
    ["2615_2617\ndelAAG", "green", "blue", "f"],
    ["1707delT", "green", "black", "f"],
    ["3183G>A", "black", "red", "b"],
    ["Duplication", "black", "red", "b"]
], columns=["Gene", "Wild", "Mutant", "MutPos"])
# ====================================================데이터 정의=====================================================


def create_moving_avg_std_cutoffs(signal, window_size, k):
    half_window = window_size // 2
    padded = np.pad(signal, (half_window, half_window), mode='reflect')
    cutoff = np.zeros_like(signal, dtype=float)
    for i in range(len(signal)):
        w = padded[i : i + window_size]; cutoff[i] = np.mean(w) + k * np.std(w)
    return cutoff


def _read_record(filepath):
    # BioPython이 FSA도 ABI 리더로 처리 가능하다는 전제 유지
    return SeqIO.read(filepath, "abi")


def _get_size_pred(record):
    return np.array(record.annotations['abif_raw']['SMap2'], dtype=float)


def _get_raw_channels(record):
    channels = {}
    for color, key in zip(COLORS, CHANNEL_KEYS):
        channels[color] = np.array(record.annotations['abif_raw'][key], dtype=float)
    return channels


def _detect_peaks_generic(filepath, best_params_dict, mode="test", return_cutoff=False):
    """
    mode: "test" | "control"
    - test: 기존 detect_peaks 동작 (global baseline 300)
    - control: 기존 detect_peaks_control 동작 (k/c/dynamic, baseline >= 500)
    """
    record = _read_record(filepath)
    size_pred = _get_size_pred(record)
    raw_channels = _get_raw_channels(record)

    # 패널 추정 (기존 startswith("s1") 보완)
    path_lower = filepath.lower()
    is_panel1 = (path_lower.startswith("s1") or "/s1/" in path_lower or "s1-" in path_lower or "panel1" in path_lower)

    # K/C factors
    K_FACTORS_P1 = [best_params_dict['P1_blue']['k_blue'], best_params_dict['P1_green']['k_green'],
                    best_params_dict['P1_black']['k_black'], best_params_dict['P1_red']['k_red']]
    C_FACTORS_P1 = [best_params_dict['P1_blue']['c_blue'], best_params_dict['P1_green']['c_green'],
                    best_params_dict['P1_black']['c_black'], best_params_dict['P1_red']['c_red']]
    K_FACTORS_P2 = [best_params_dict['P2_blue']['k_blue'], best_params_dict['P2_green']['k_green'],
                    best_params_dict['P2_black']['k_black'], best_params_dict['P2_red']['k_red']]
    C_FACTORS_P2 = [best_params_dict['P2_blue']['c_blue'], best_params_dict['P2_green']['c_green'],
                    best_params_dict['P2_black']['c_black'], best_params_dict['P2_red']['c_red']]

    all_peaks, cutoff_lines = [], {}

    for ind, (color, raw_signal) in enumerate(raw_channels.items()):
        smoothed_signal = uniform_filter1d(raw_signal, 7)
        global_mean, global_std = np.mean(smoothed_signal), np.std(smoothed_signal)

        if mode == "control":
            if is_panel1:
                k, c = K_FACTORS_P1[ind], C_FACTORS_P1[ind]
            else:
                k, c = K_FACTORS_P2[ind], C_FACTORS_P2[ind]

            global_baseline = max(c + global_mean + k * global_std, 500)
            dynamic_line = create_moving_avg_std_cutoffs(smoothed_signal, WINDOW_SIZE, k)
            final_cutoff_line = np.maximum(global_baseline, dynamic_line)
            peaks, _ = find_peaks(smoothed_signal, height=final_cutoff_line, distance=20, prominence=100)
        else:
            # test 모드: 기존 로직 유지
            final_cutoff_line = np.full_like(smoothed_signal, max(global_mean + 1.5 * global_std, 300))
            peaks, _ = find_peaks(smoothed_signal, height=final_cutoff_line)

        for idx in peaks:
            all_peaks.append((color, size_pred[idx], smoothed_signal[idx]))
        cutoff_lines[color] = (size_pred, final_cutoff_line)

    if return_cutoff:
        return all_peaks, cutoff_lines
    return all_peaks


def detect_peaks(filepath, best_params_dict, return_cutoff=False):
    return _detect_peaks_generic(filepath, best_params_dict, mode="test", return_cutoff=return_cutoff)

def detect_peaks_control(filepath, best_params_dict, return_cutoff=False):
    return _detect_peaks_generic(filepath, best_params_dict, mode="control", return_cutoff=return_cutoff)


def get_windows(file_path, best_params, mutation_map, window_size=2, iscontrol = False):
    if iscontrol:
        peaks = detect_peaks_control(file_path, best_params)
    else:
        # 오타 수정: best_parms -> best_params
        peaks = detect_peaks(file_path, best_params)

    windows = []
    for idx, (_, bp_pos, _) in enumerate(peaks):
        if idx >= len(mutation_map):
            break

        mutpos = mutation_map.iloc[idx]["MutPos"]

        start = bp_pos - (window_size + 0.7)
        end = bp_pos + (window_size + 0.7)

        if mutpos == "f":
            start -= 0.7
        elif mutpos == "b":
            end += 0.7

        windows.append((start, end))

    return windows



def fix_windows(windows, control_peaks=None, mutation_map=None):
    """
    Modify windows to connect front and back,
    adjusting the boundary position based on the ‘f’/'b' combination.
    - f (front): Shorten the front window (ratio↓)
    - b (back): Lengthen the front window (ratio↑)
    - f-b or b-f combination: Midpoint value (0.5)
    """
    windows = sorted(windows, key=lambda x: x[0])

    if control_peaks is None or mutation_map is None:
        raise ValueError("control_peaks와 mutation_map이 모두 필요합니다.")

    control_peaks = sorted(control_peaks, key=lambda x: x[1])
    fixed = []

    for i in range(len(control_peaks)):
        curr_peak = control_peaks[i][1]

        if i == 0:
            start = windows[0][0]
        else:
            start = fixed[-1][1]

        # 마지막 윈도우 예외
        if i >= len(control_peaks) - 1:
            end = curr_peak + 5
            fixed.append((start, end))
            continue

        next_peak = control_peaks[i + 1][1]
        curr_pos = mutation_map.iloc[i]["MutPos"]
        next_pos = mutation_map.iloc[i + 1]["MutPos"]

        # --- 비율 결정 ---
        if curr_pos == "b" and next_pos == "f":
            ratio = 0.5        # 양쪽 타협
        elif curr_pos == "f" and next_pos == "b":
            ratio = 0.5        # 반대 방향도 동일하게 타협
        elif curr_pos == "b" and next_pos == "b":
            ratio = 0.70       # 둘 다 뒤쪽 강조
        elif curr_pos == "f" and next_pos == "f":
            ratio = 0.30       # 둘 다 앞쪽 강조
        elif curr_pos == "b":
            ratio = 0.7
        elif curr_pos == "f":
            ratio = 0.3
        else:
            ratio = 0.5

        # --- 경계 계산 ---
        end = curr_peak + (next_peak - curr_peak) * ratio
        fixed.append((start, end))

    return fixed



def get_control_peak_positions(control_path, best_params, control_windows, mutation_map):
    """
    컨트롤 패널에서 각 윈도우의 대표 피크 위치를 추출
    
    Returns:
        List of dict: [{'wild_pos': float, 'mutant_pos': float, 'colors': [wild, mutant]}, ...]
    """
    peaks = detect_peaks_control(control_path, best_params)
    peaks = sorted(peaks, key=lambda x: x[1])
    
    control_positions = []
    
    for idx, (s, e) in enumerate(control_windows):
        in_window = [(c, p, inten) for c, p, inten in peaks if s <= p <= e]
        wild, mut, pos = mutation_map.iloc[idx][["Wild", "Mutant", "MutPos"]]
        
        # 허용된 색만 고려
        filtered = [x for x in in_window if x[0] in [wild, mut]]
        
        if not filtered:
            control_positions.append({
                'wild_pos': None, 
                'mutant_pos': None, 
                'colors': [wild, mut],
                'mutpos': pos
            })
            continue
        
        # 위치 기준 정렬
        filtered = sorted(filtered, key=lambda x: x[1])
        
        # wild와 mutant 위치 추출
        wild_peaks = [x for x in filtered if x[0] == wild]
        mut_peaks = [x for x in filtered if x[0] == mut]
        
        wild_pos = wild_peaks[0][1] if wild_peaks else None
        mut_pos = mut_peaks[0][1] if mut_peaks else None
        
        control_positions.append({
            'wild_pos': wild_pos,
            'mutant_pos': mut_pos,
            'colors': [wild, mut],
            'mutpos': pos
        })
    
    return control_positions


def show_with_windows(file_path, best_params, windows=None, title="Test panel", mutation_map=None, visual=False):
    record = _read_record(file_path)
    size_pred = _get_size_pred(record)

    if visual:
        plt.figure(figsize=(20, 8))

        trace = {}
        for color, key in zip(COLORS, CHANNEL_KEYS):
            raw = np.array(record.annotations['abif_raw'][key], dtype=float)
            trace[color] = uniform_filter1d(raw, 7)
            plt.plot(size_pred, trace[color], color=color, label=f"{color}")

        peaks = detect_peaks(file_path, best_params)
        for (color, bp_pos, intensity) in peaks:
            plt.scatter(bp_pos, intensity, color="black", s=40, marker="x")

        if windows is not None:
            for idx, (start, end) in enumerate(windows):
                plt.axvspan(start, end, color="gray", alpha=0.2)
                if mutation_map is not None and idx < len(mutation_map):
                    gene = mutation_map.iloc[idx]["Gene"]
                    center = (start + end) / 2
                    plt.text(center, plt.ylim()[1]*1.05, gene, 
                            ha="center", va="top", fontsize=10, rotation=0, color="black")

        plt.title(title)
        plt.xlim([20, 80])
        ymin, ymax = plt.ylim()
        plt.ylim(ymin, ymax * 1.2)
        plt.legend()
        plt.show()

        return peaks
    else:
        # 시각화가 아닐 때는 피크만 반환
        return detect_peaks(file_path, best_params)


def evaluate_panel(test_path, control_windows, mutation_map, control_positions=None, 
                   position_tolerance=0.3, visual_f=False, ratio = 0.4):
    """
    Args:
        position_tolerance: 컨트롤 피크 위치로부터 허용 오차 (bp 단위)
    """
    peaks = show_with_windows(test_path, all_best_params,
                              windows=control_windows,
                              title="Test with Control Windows", 
                              mutation_map=mutation_map, visual=visual_f)

    peaks = sorted(peaks, key=lambda x: x[1])

    results = []
    for idx, (s, e) in enumerate(control_windows):
        in_window = [(c, inten, p) for c, p, inten in peaks if s <= p <= e]

        wild, mut, pos = mutation_map.iloc[idx][["Wild", "Mutant", "MutPos"]]

        # 1) 허용된 색만 고려
        filtered = [x for x in in_window if x[0] in [wild, mut]]

        # 2) 🔹 컨트롤 피크 위치 기반 필터링 추가
        if control_positions is not None and idx < len(control_positions):
            ctrl_info = control_positions[idx]
            wild_ref = ctrl_info['wild_pos']
            mut_ref = ctrl_info['mutant_pos']
            
            position_filtered = []
            for color, inten, peak_pos in filtered:
                # 해당 색상의 기준 위치 선택
                ref_pos = wild_ref if color == wild else mut_ref
                
                # 기준 위치가 없으면 통과
                if ref_pos is None:
                    position_filtered.append((color, inten, peak_pos))
                    continue
                
                # 위치 차이 계산
                pos_diff = abs(peak_pos - ref_pos)
                
                # tolerance 내에 있으면 유지
                if pos_diff <= position_tolerance:
                    position_filtered.append((color, inten, peak_pos))
            
            filtered = position_filtered

        # 3) 강도 기반 필터링
        if filtered:
            # 최고 intensity 구하기
            max_inten = max(inten for _, inten, _ in filtered)

            # ratio 기준으로 우선 필터링
            filtered = [x for x in filtered if x[1] >= ratio * max_inten]

            # 🔹 최고 intensity 피크를 별도로 보존
            top_peak = max(filtered, key=lambda x: x[1])

            # 🔹 나머지 피크 중 intensity < 500은 제외
            filtered = [
                x for x in filtered
                if (x == top_peak) or (x[1] >= 400)
            ]

        filtered = sorted(filtered, key=lambda x: x[2])
        colors_present = [c for c, _, _ in filtered]

        # 4) 판정 로직
        if len(colors_present) == 0:
            call = 0
        elif len(colors_present) == 1:
            call = 0 if colors_present[0] == wild else 2
        elif len(colors_present) >= 2:
            if pos == "b" and colors_present[:2] == [wild, mut]:
                call = 1
            elif pos == "f" and colors_present[:2] == [mut, wild]:
                call = 1
            else:
                call = 0
        else:
            call = -1

        results.append(call)

    return results


def check_diplotype(a1_name, a2_name, sample, ref, alleles_ref):
    """diplotype 조합이 sample과 일치하는지 확인"""
    if a1_name not in alleles_ref or a2_name not in alleles_ref:
        return False
    
    a1_seq = alleles_ref[a1_name]
    a2_seq = alleles_ref[a2_name]
    
    if a1_seq is None or a2_seq is None:
        return False
    
    has_star5 = (a1_name == "*5" or a2_name == "*5")
    is_star5_homo = (a1_name == "*5" and a2_name == "*5")
    
    for position, genotype_code in enumerate(sample):
        ref_base = ref[position]
        a1_base = a1_seq[position]
        a2_base = a2_seq[position]
        
        if has_star5:
            if is_star5_homo:
                if position == 9:
                    if genotype_code != 2:
                        return False
                else:
                    if genotype_code != 0:
                        return False
            else:
                other_seq = a2_seq if a1_name == "*5" else a1_seq
                other_base = other_seq[position]
                
                if position == 9:
                    if other_base == ref_base:
                        if genotype_code != 1:
                            return False
                    else:
                        if genotype_code != 1:
                            return False
                else:
                    if other_base == ref_base:
                        if genotype_code != 0:
                            return False
                    else:
                        if genotype_code != 2:
                            return False
            continue
        
        if a1_base == '-' or a2_base == '-':
            return False
        
        expected_code = calculate_genotype_code(a1_base, a2_base, ref_base)
        if genotype_code != expected_code:
            return False
    
    return True


def calculate_genotype_code(a1_base, a2_base, ref_base):
    """두 allele base로부터 genotype code 계산"""
    a1_is_ref = (a1_base == ref_base)
    a2_is_ref = (a2_base == ref_base)
    
    if a1_is_ref and a2_is_ref:
        return 0
    elif a1_is_ref or a2_is_ref:
        return 1
    elif a1_base == a2_base:
        return 2
    else:
        return 1


def sort_allele_name(allele_name):
    """Allele 이름을 정렬 가능한 키로 변환"""
    if '(or ' in allele_name:
        allele_name = allele_name.split('(or ')[0].strip()
    
    name = allele_name.strip('*')
    
    num_part = ''
    suffix = ''
    for i, char in enumerate(name):
        if char.isdigit():
            num_part += char
        else:
            suffix = name[i:]
            break
    
    num = int(num_part) if num_part else float('inf')
    
    return (num, suffix)


def call_diplotype(sample, alleles_ref):
    ref = alleles_ref["*1"]
    candidates = []

    has_10D = (sample[1] in (1,2) and sample[9] in (1,2))
    has_5 = (sample[1] not in (1, 2) and sample[9] in (1, 2))

    if has_10D:
        reduced_sample = [v if i not in [1,9] else 0 for i, v in enumerate(sample)]
        
        for idx, i in enumerate(reduced_sample):
            if i == 2:
                reduced_sample[idx] = 1
                
        for allele in alleles_ref.keys():
            if allele in ["*5", "*10B", "*10D"]:
                continue
            if check_diplotype(allele, "*1", reduced_sample, ref, alleles_ref):
                sorted_pair = sorted([allele, "*10D"], key=sort_allele_name)
                candidates.append(f"{sorted_pair[0]}/{sorted_pair[1]}")
                if allele == "*1":
                    candidates.append(f"*5/*10B")

    if has_5:
        reduced_sample = [v if i not in [9] else 0 for i, v in enumerate(sample)]
        
        for idx, i in enumerate(reduced_sample):
            if i == 2:
                reduced_sample[idx] = 1
                
        for allele in alleles_ref.keys():
            if allele in ["*5", "*10B", "*10D"]:
                continue
            if check_diplotype(allele, "*1", reduced_sample, ref, alleles_ref):
                sorted_pair = sorted([allele, "*5"], key=sort_allele_name)
                candidates.append(f"{sorted_pair[0]}/{sorted_pair[1]}")

    if not has_5 and not has_10D:
        for a1, a2 in itertools.combinations_with_replacement(alleles_ref.keys(), 2):
            if check_diplotype(a1, a2, sample, ref, alleles_ref):
                sorted_pair = sorted([a1, a2], key=sort_allele_name)
                candidates.append(f"{sorted_pair[0]}/{sorted_pair[1]}")

    candidates = sorted(list(set(candidates)),
                        key=lambda x: (sort_allele_name(x.split('/')[0]), sort_allele_name(x.split('/')[1])))
    return candidates
