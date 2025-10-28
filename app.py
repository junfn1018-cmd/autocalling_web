from flask import Flask, render_template, request, jsonify, send_file
import os
import tempfile
import json
import logging
from werkzeug.utils import secure_filename
import matplotlib
matplotlib.use('Agg')  # GUI 없이 사용
import matplotlib.pyplot as plt
import io
import base64

# 로깅 설정
logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)
try:
    from autocalling_web import (
        get_windows, fix_windows, get_control_peak_positions, 
        evaluate_panel, call_diplotype, all_best_params, 
        alleles_ref, mutation_map
    )
    print("autocalling_web 모듈 import 성공")
except Exception as e:
    print(f"autocalling_web 모듈 import 실패: {e}")
    import traceback
    traceback.print_exc()

app = Flask(__name__)
app.config['MAX_CONTENT_LENGTH'] = 16 * 1024 * 1024  # 16MB max file size
app.config['UPLOAD_FOLDER'] = 'uploads'

# 업로드 폴더 생성
os.makedirs(app.config['UPLOAD_FOLDER'], exist_ok=True)

ALLOWED_EXTENSIONS = {'abi', 'fsa'}

def allowed_file(filename):
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/test')
def test_endpoint():
    print("=== 테스트 엔드포인트 호출됨 ===")
    return jsonify({'status': 'success', 'message': '테스트 성공'})

@app.route('/simple-upload', methods=['POST'])
def simple_upload():
    print("=== 간단한 업로드 테스트 ===")
    try:
        print(f"Request method: {request.method}")
        print(f"Request files: {list(request.files.keys())}")
        return jsonify({'status': 'success', 'message': '업로드 테스트 성공'})
    except Exception as e:
        print(f"간단한 업로드 오류: {e}")
        return jsonify({'error': str(e)}), 500

@app.route('/upload', methods=['POST'])
def upload_files():
    logger.info("=== 함수 진입 ===")
    logger.info(f"Request method: {request.method}")
    logger.info(f"Request files: {list(request.files.keys())}")
    try:
        logger.info("=== 업로드 요청 시작 ===")
        
        # 파일 업로드 확인
        if ('control_panel1' not in request.files or 'control_panel2' not in request.files or 
            'test_panel1' not in request.files or 'test_panel2' not in request.files):
            logger.info("오류: 모든 패널 파일이 없음")
            return jsonify({'error': '모든 패널 파일이 필요합니다.'}), 400
        
        control_panel1 = request.files['control_panel1']
        control_panel2 = request.files['control_panel2']
        test_panel1 = request.files['test_panel1']
        test_panel2 = request.files['test_panel2']
        
        logger.info(f"컨트롤 패널1: {control_panel1.filename}")
        logger.info(f"컨트롤 패널2: {control_panel2.filename}")
        logger.info(f"테스트 패널1: {test_panel1.filename}")
        logger.info(f"테스트 패널2: {test_panel2.filename}")
        
        if (control_panel1.filename == '' or control_panel2.filename == '' or 
            test_panel1.filename == '' or test_panel2.filename == ''):
            logger.info("오류: 파일명이 비어있음")
            return jsonify({'error': '모든 패널 파일을 선택해주세요.'}), 400
        
        if not (allowed_file(control_panel1.filename) and allowed_file(control_panel2.filename) and
                allowed_file(test_panel1.filename) and allowed_file(test_panel2.filename)):
            logger.info("오류: 지원하지 않는 파일 형식")
            return jsonify({'error': 'ABI 또는 FSA 파일만 업로드 가능합니다.'}), 400
        
        # 파일 저장
        control_panel1_filename = secure_filename(control_panel1.filename)
        control_panel2_filename = secure_filename(control_panel2.filename)
        test_panel1_filename = secure_filename(test_panel1.filename)
        test_panel2_filename = secure_filename(test_panel2.filename)
        
        control_panel1_path = os.path.join(app.config['UPLOAD_FOLDER'], control_panel1_filename)
        control_panel2_path = os.path.join(app.config['UPLOAD_FOLDER'], control_panel2_filename)
        test_panel1_path = os.path.join(app.config['UPLOAD_FOLDER'], test_panel1_filename)
        test_panel2_path = os.path.join(app.config['UPLOAD_FOLDER'], test_panel2_filename)
        
        # 디버그: 파일 경로 출력
        logger.info(f"=== 파일 업로드 디버그 ===")
        logger.info(f"컨트롤 패널1 경로: {control_panel1_path}")
        logger.info(f"컨트롤 패널2 경로: {control_panel2_path}")
        logger.info(f"테스트 패널1 경로: {test_panel1_path}")
        logger.info(f"테스트 패널2 경로: {test_panel2_path}")
        
        control_panel1.save(control_panel1_path)
        control_panel2.save(control_panel2_path)
        test_panel1.save(test_panel1_path)
        test_panel2.save(test_panel2_path)
        
        # 디버그: 파일 저장 후 확인
        logger.info(f"=== 파일 저장 후 확인 ===")
        logger.info(f"컨트롤 패널1 저장됨: {os.path.exists(control_panel1_path)}")
        logger.info(f"컨트롤 패널2 저장됨: {os.path.exists(control_panel2_path)}")
        logger.info(f"테스트 패널1 저장됨: {os.path.exists(test_panel1_path)}")
        logger.info(f"테스트 패널2 저장됨: {os.path.exists(test_panel2_path)}")
        
        # 분석 수행
        logger.info("analyze_samples 함수 호출 시작")
        result = analyze_samples(control_panel1_path, control_panel2_path, test_panel1_path, test_panel2_path)
        logger.info("analyze_samples 함수 호출 완료")
        
        # 임시 파일 삭제
        os.remove(control_panel1_path)
        os.remove(control_panel2_path)
        os.remove(test_panel1_path)
        os.remove(test_panel2_path)
        
        return jsonify(result)
        
    except Exception as e:
        import traceback
        error_details = traceback.format_exc()
        print(f"상세 오류 정보:\n{error_details}")
        return jsonify({
            'error': f'분석 중 오류가 발생했습니다: {str(e)}',
            'details': error_details
        }), 500

def create_visualization(file_path, windows, mutation_map, title="Analysis Result"):
    """분석 결과를 시각화하여 base64 이미지로 반환"""
    try:
        from Bio import SeqIO
        import numpy as np
        from scipy.ndimage import uniform_filter1d
        
        record = SeqIO.read(file_path, "abi")
        size_pred = np.array(record.annotations['abif_raw']['SMap2'], dtype=float)
        
        colors = ['blue', 'green', 'black', 'red']
        channels = ['DATA9', 'DATA10', 'DATA11', 'DATA12']
        
        plt.figure(figsize=(20, 8))
        
        # 각 채널의 데이터 플롯
        for channel, color in zip(channels, colors):
            raw = np.array(record.annotations['abif_raw'][channel], dtype=float)
            trace = uniform_filter1d(raw, 7)
            plt.plot(size_pred, trace, color=color, label=f"{color}")
        
        # 윈도우 표시
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
        plt.xlabel('Size (bp)')
        plt.ylabel('Intensity')
        
        # 이미지를 base64로 변환
        img_buffer = io.BytesIO()
        plt.savefig(img_buffer, format='png', dpi=100, bbox_inches='tight')
        img_buffer.seek(0)
        img_base64 = base64.b64encode(img_buffer.getvalue()).decode()
        plt.close()
        
        return img_base64
        
    except Exception as e:
        logger.error(f"시각화 생성 오류: {e}")
        return None

def analyze_samples(control_panel1_path, control_panel2_path, test_panel1_path, test_panel2_path):
    """샘플 분석을 수행하는 함수 - 4개 패널 파일 처리"""
    logger.info("analyze_samples 함수 시작")
    try:
        logger.info(f"분석 시작 - 컨트롤 패널1: {control_panel1_path}")
        logger.info(f"분석 시작 - 컨트롤 패널2: {control_panel2_path}")
        logger.info(f"분석 시작 - 테스트 패널1: {test_panel1_path}")
        logger.info(f"분석 시작 - 테스트 패널2: {test_panel2_path}")
        
        # 파일 존재 확인
        files = [control_panel1_path, control_panel2_path, test_panel1_path, test_panel2_path]
        for file_path in files:
            if not os.path.exists(file_path):
                raise FileNotFoundError(f"파일을 찾을 수 없습니다: {file_path}")
        
        logger.info("1단계: 컨트롤 윈도우 생성 중...")
        # 1. 컨트롤 윈도우 생성 (패널1과 패널2 모두 사용)
        control_windows1 = get_windows(control_panel1_path, all_best_params, mutation_map, iscontrol=True)
        control_windows2 = get_windows(control_panel2_path, all_best_params, mutation_map, iscontrol=True)
        logger.info(f"컨트롤 윈도우 생성 완료 - 패널1: {len(control_windows1)}개, 패널2: {len(control_windows2)}개")
        
        logger.info("2단계: mutation_map 분할 중...")
        # 원본 코드 로직: mutation_map을 패널별로 분할
        mutation_map_p1 = mutation_map.iloc[:10].reset_index(drop=True)  # 패널1: 0-9
        mutation_map_p2 = mutation_map.iloc[10:].reset_index(drop=True)   # 패널2: 10-17
        logger.info(f"mutation_map 분할 완료 - 패널1: {len(mutation_map_p1)}개, 패널2: {len(mutation_map_p2)}개")
        
        logger.info("3단계: 윈도우 길이 조정 중...")
        # 윈도우 길이를 각 패널의 mutation_map 길이에 맞춤
        control_windows1 = control_windows1[:len(mutation_map_p1)]
        control_windows2 = control_windows2[:len(mutation_map_p2)]
        logger.info(f"윈도우 길이 조정 완료 - 패널1: {len(control_windows1)}개, 패널2: {len(control_windows2)}개")
        
        logger.info("4단계: 컨트롤 피크 추출 중...")
        # 원본 코드 로직: detect_peaks_control 사용
        from autocalling_web import detect_peaks_control
        control_peaks1 = detect_peaks_control(control_panel1_path, all_best_params)
        control_peaks2 = detect_peaks_control(control_panel2_path, all_best_params)
        logger.info(f"컨트롤 피크 추출 완료 - 패널1: {len(control_peaks1)}개, 패널2: {len(control_peaks2)}개")
        
        logger.info("5단계: 윈도우 수정 중...")
        # 원본 코드 로직: fix_windows 사용
        fixed_windows1 = fix_windows(control_windows1, control_peaks1, mutation_map_p1)
        fixed_windows2 = fix_windows(control_windows2, control_peaks2, mutation_map_p2)
        logger.info(f"윈도우 수정 완료 - 패널1: {len(fixed_windows1)}개, 패널2: {len(fixed_windows2)}개")
        
        logger.info("6단계: 컨트롤 피크 위치 추출 중...")
        # 컨트롤 피크 위치 추출
        control_positions1 = get_control_peak_positions(control_panel1_path, all_best_params, fixed_windows1, mutation_map_p1)
        control_positions2 = get_control_peak_positions(control_panel2_path, all_best_params, fixed_windows2, mutation_map_p2)
        logger.info(f"컨트롤 피크 위치 추출 완료 - 패널1: {len(control_positions1)}개, 패널2: {len(control_positions2)}개")
        
        logger.info("7단계: 테스트 샘플 분석 중...")
        # 원본 코드 로직: ratio를 조절하면서 여러 번 분석
        ratios = [0.50, 0.39, 0.28]  # 원본 코드와 동일한 ratio 시퀀스
        all_diplotype_results = []
        
        for ratio_idx, ratio in enumerate(ratios):
            logger.info(f"분석 {ratio_idx + 1}/3 - ratio: {ratio}")
            
            # 각 패널의 mutation_map과 윈도우 사용
            test_results1 = evaluate_panel(
                test_panel1_path, 
                fixed_windows1, 
                mutation_map_p1, 
                control_positions=control_positions1,
                visual_f=False,
                ratio=ratio
            )
            test_results2 = evaluate_panel(
                test_panel2_path, 
                fixed_windows2, 
                mutation_map_p2, 
                control_positions=control_positions2,
                visual_f=False,
                ratio=ratio
            )
            
            # 원본 코드: sample = p1 + p2
            combined_results = test_results1 + test_results2
            
            # 다이플로타입 호출
            diplotype_candidates = call_diplotype(combined_results, alleles_ref)
            
            logger.info(f"ratio {ratio} 결과: {diplotype_candidates}")
            all_diplotype_results.append({
                'ratio': ratio,
                'diplotypes': diplotype_candidates,
                'test_results': combined_results
            })
        
        logger.info(f"모든 분석 완료 - 총 {len(all_diplotype_results)}번의 분석")
        
        logger.info("8단계: 결과 비교 및 정리 중...")
        # 모든 결과가 동일한지 확인
        first_result = all_diplotype_results[0]['diplotypes']
        all_same = all(result['diplotypes'] == first_result for result in all_diplotype_results)
        
        logger.info(f"첫 번째 결과: {first_result}")
        logger.info(f"모든 결과 동일 여부: {all_same}")
        
        if all_same:
            # 모든 결과가 동일하면 하나의 결과로 반환
            logger.info("모든 ratio에서 동일한 결과 - 단일 결과 반환")
            final_diplotype_candidates = first_result
            final_test_results = all_diplotype_results[0]['test_results']
            ratio_analysis = None
        else:
            # 결과가 다르면 각 ratio별 결과를 제시
            logger.info("ratio별로 다른 결과 - 각 ratio별 결과 제시")
            # 첫 번째 결과를 기본으로 사용하되, 빈 리스트가 아닌 경우에만
            final_diplotype_candidates = first_result if first_result else []
            final_test_results = all_diplotype_results[0]['test_results']
            ratio_analysis = []
            
            for result in all_diplotype_results:
                ratio_analysis.append({
                    'ratio': result['ratio'],
                    'diplotypes': result['diplotypes'],
                    'count': len(result['diplotypes'])
                })
                logger.info(f"ratio {result['ratio']}: {result['diplotypes']}")
        
        logger.info(f"최종 다이플로타입 후보: {final_diplotype_candidates}")
        logger.info(f"최종 다이플로타입 후보 개수: {len(final_diplotype_candidates)}")
        
        # 9. 시각화 생성
        logger.info("시각화 생성 중...")
        visualization_panel1 = create_visualization(
            test_panel1_path, fixed_windows1, mutation_map_p1, 
            title="Test Panel 1 Analysis"
        )
        visualization_panel2 = create_visualization(
            test_panel2_path, fixed_windows2, mutation_map_p2, 
            title="Test Panel 2 Analysis"
        )
        logger.info("시각화 생성 완료")
        
        # 10. 결과 정리
        mutation_names = mutation_map['Gene'].tolist()
        
        analysis_result = {
            'success': True,
            'test_results': final_test_results,
            'test_results_panel1': final_test_results[:10],  # 패널1 결과
            'test_results_panel2': final_test_results[10:],   # 패널2 결과
            'mutation_names': mutation_names,
            'diplotype_candidates': final_diplotype_candidates,
            'ratio_analysis': ratio_analysis,
            'visualization_panel1': visualization_panel1,
            'visualization_panel2': visualization_panel2,
            'summary': {
                'total_mutations': len(final_test_results),
                'detected_mutations': len([r for r in final_test_results if r > 0]),
                'diplotype_count': len(final_diplotype_candidates),
                'analysis_count': len(all_diplotype_results),
                'consistent_results': all_same
            }
        }
        
        logger.info("분석 완료!")
        return analysis_result
        
    except Exception as e:
        import traceback
        error_details = traceback.format_exc()
        print(f"analyze_samples 오류:\n{error_details}")
        raise Exception(f"분석 실패: {str(e)}")

@app.route('/api/alleles')
def get_alleles():
    """사용 가능한 대립유전자 목록 반환"""
    return jsonify(list(alleles_ref.keys()))

@app.route('/api/mutations')
def get_mutations():
    """돌연변이 맵 정보 반환"""
    mutations = mutation_map.to_dict('records')
    return jsonify(mutations)

if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0', port=5004)
