// DOM 요소들
const uploadForm = document.getElementById('uploadForm');
const analyzeBtn = document.getElementById('analyzeBtn');
const loadingSpinner = document.getElementById('loadingSpinner');
const resultsSection = document.getElementById('resultsSection');
const errorAlert = document.getElementById('errorAlert');
const errorMessage = document.getElementById('errorMessage');

// 결과 표시 요소들
const totalMutations = document.getElementById('totalMutations');
const detectedMutations = document.getElementById('detectedMutations');
const diplotypeCount = document.getElementById('diplotypeCount');
const diplotypeResults = document.getElementById('diplotypeResults');
const resultsTable = document.getElementById('resultsTable');
const visualizationPanel1 = document.getElementById('visualizationPanel1');
const visualizationPanel2 = document.getElementById('visualizationPanel2');

// 폼 제출 이벤트 리스너
uploadForm.addEventListener('submit', async function(e) {
    e.preventDefault();
    
    // 파일 유효성 검사
    const controlPanel1 = document.getElementById('controlPanel1').files[0];
    const controlPanel2 = document.getElementById('controlPanel2').files[0];
    const testPanel1 = document.getElementById('testPanel1').files[0];
    const testPanel2 = document.getElementById('testPanel2').files[0];
    
    if (!controlPanel1 || !controlPanel2 || !testPanel1 || !testPanel2) {
        showError('모든 패널 파일을 선택해주세요.');
        return;
    }
    
    if (!isValidABIFile(controlPanel1) || !isValidABIFile(controlPanel2) || 
        !isValidABIFile(testPanel1) || !isValidABIFile(testPanel2)) {
        showError('ABI 또는 FSA 파일만 업로드 가능합니다.');
        return;
    }
    
    // 분석 시작
    await startAnalysis();
});

// ABI 파일 유효성 검사
function isValidABIFile(file) {
    const allowedExtensions = ['.abi', '.fsa'];
    const fileName = file.name.toLowerCase();
    return allowedExtensions.some(ext => fileName.endsWith(ext));
}

// 분석 시작
async function startAnalysis() {
    try {
        // UI 상태 변경
        setLoadingState(true);
        hideError();
        hideResults();
        
        // FormData 생성
        const formData = new FormData();
        formData.append('control_panel1', document.getElementById('controlPanel1').files[0]);
        formData.append('control_panel2', document.getElementById('controlPanel2').files[0]);
        formData.append('test_panel1', document.getElementById('testPanel1').files[0]);
        formData.append('test_panel2', document.getElementById('testPanel2').files[0]);
        
        // API 호출 - 먼저 간단한 테스트
        console.log('간단한 업로드 테스트 시작');
        const testResponse = await fetch('/simple-upload', {
            method: 'POST',
            body: formData
        });
        console.log('간단한 테스트 응답:', testResponse.status);
        
        if (!testResponse.ok) {
            throw new Error('간단한 업로드 테스트 실패');
        }
        
        console.log('API 호출 시작');
        const response = await fetch('/upload', {
            method: 'POST',
            body: formData
        });
        console.log('API 응답 받음:', response.status);
        
        const result = await response.json();
        
        if (!response.ok) {
            throw new Error(result.error || '분석 중 오류가 발생했습니다.');
        }
        
        if (result.success) {
            displayResults(result);
        } else {
            throw new Error('분석이 실패했습니다.');
        }
        
    } catch (error) {
        console.error('Analysis error:', error);
        showError(error.message);
    } finally {
        setLoadingState(false);
    }
}

// 로딩 상태 설정
function setLoadingState(loading) {
    if (loading) {
        analyzeBtn.disabled = true;
        analyzeBtn.innerHTML = '<i class="fas fa-spinner fa-spin me-2"></i>분석 중...';
        loadingSpinner.style.display = 'block';
    } else {
        analyzeBtn.disabled = false;
        analyzeBtn.innerHTML = '<i class="fas fa-play me-2"></i>분석 시작';
        loadingSpinner.style.display = 'none';
    }
}

// 오류 메시지 표시
function showError(message) {
    errorMessage.textContent = message;
    errorAlert.style.display = 'block';
    errorAlert.scrollIntoView({ behavior: 'smooth' });
}

// 오류 메시지 숨기기
function hideError() {
    errorAlert.style.display = 'none';
}

// 결과 숨기기
function hideResults() {
    resultsSection.style.display = 'none';
}

// 결과 표시
function displayResults(result) {
    // 요약 정보 업데이트
    updateSummary(result.summary);
    
    // 다이플로타입 결과 표시
    displayDiplotypeResults(result.diplotype_candidates, result.ratio_analysis);
    
    // 상세 결과 테이블 표시
    displayDetailedResults(result.test_results, result.mutation_names);
    
    // 시각화 결과 표시
    displayVisualizations(result.visualization_panel1, result.visualization_panel2);
    
    // 결과 섹션 표시
    resultsSection.style.display = 'block';
    resultsSection.classList.add('fade-in');
    
    // 결과 섹션으로 스크롤
    resultsSection.scrollIntoView({ behavior: 'smooth' });
}

// 요약 정보 업데이트
function updateSummary(summary) {
    totalMutations.textContent = summary.total_mutations;
    detectedMutations.textContent = summary.detected_mutations;
    diplotypeCount.textContent = summary.diplotype_count;
}

// 다이플로타입 결과 표시
function displayDiplotypeResults(diplotypeCandidates, ratioAnalysis) {
    diplotypeResults.innerHTML = '';
    
    if (diplotypeCandidates.length === 0) {
        diplotypeResults.innerHTML = `
            <div class="text-center text-muted py-4">
                <i class="fas fa-exclamation-circle fa-2x mb-2"></i>
                <p>검출된 다이플로타입이 없습니다.</p>
            </div>
        `;
        return;
    }
    
    // ratio 분석 결과가 있는 경우 (결과가 다른 경우)
    if (ratioAnalysis && ratioAnalysis.length > 0) {
        diplotypeResults.innerHTML = `
            <div class="alert alert-warning mb-3">
                <i class="fas fa-info-circle me-2"></i>
                <strong>주의:</strong> 서로 다른 ratio에서 다른 결과가 나타났습니다.
            </div>
        `;
        
        // 각 ratio별 결과 표시
        ratioAnalysis.forEach((analysis, index) => {
            const ratioCard = document.createElement('div');
            ratioCard.className = 'card mb-3';
            ratioCard.innerHTML = `
                <div class="card-header bg-light">
                    <h6 class="mb-0">
                        <i class="fas fa-chart-line me-2"></i>
                        Ratio ${analysis.ratio} 분석 결과
                    </h6>
                </div>
                <div class="card-body">
                    <div class="row">
                        <div class="col-md-8">
                            ${analysis.diplotypes.map(diplotype => `
                                <div class="diplotype-item slide-up mb-2" style="animation-delay: ${index * 0.1}s">
                                    <h6><i class="fas fa-dna me-2"></i>${diplotype}</h6>
                                    <p class="mb-0 text-muted">가능한 유전자형 조합</p>
                                </div>
                            `).join('')}
                        </div>
                        <div class="col-md-4 text-center">
                            <div class="badge bg-primary fs-6">${analysis.count}개 후보</div>
                        </div>
                    </div>
                </div>
            `;
            diplotypeResults.appendChild(ratioCard);
        });
    } else {
        // 모든 ratio에서 동일한 결과인 경우
        diplotypeResults.innerHTML = `
            <div class="alert alert-success mb-3">
                <i class="fas fa-check-circle me-2"></i>
                <strong>일관된 결과:</strong> 모든 ratio에서 동일한 결과가 나타났습니다.
            </div>
        `;
        
        diplotypeCandidates.forEach((diplotype, index) => {
            const diplotypeItem = document.createElement('div');
            diplotypeItem.className = 'diplotype-item slide-up';
            diplotypeItem.style.animationDelay = `${index * 0.1}s`;
            
            diplotypeItem.innerHTML = `
                <h6><i class="fas fa-dna me-2"></i>${diplotype}</h6>
                <p class="mb-0">가능한 유전자형 조합</p>
            `;
            
            diplotypeResults.appendChild(diplotypeItem);
        });
    }
}

// 상세 결과 테이블 표시
function displayDetailedResults(testResults, mutationNames) {
    resultsTable.innerHTML = '';
    
    testResults.forEach((result, index) => {
        const row = document.createElement('tr');
        row.className = 'slide-up';
        row.style.animationDelay = `${index * 0.05}s`;
        
        const resultInfo = getResultInfo(result);
        
        row.innerHTML = `
            <td>${index + 1}</td>
            <td><strong>${mutationNames[index]}</strong></td>
            <td>
                <span class="result-badge ${resultInfo.class}">
                    ${resultInfo.text}
                </span>
            </td>
            <td class="interpretation">${resultInfo.interpretation}</td>
        `;
        
        resultsTable.appendChild(row);
    });
}

// 결과 정보 가져오기
function getResultInfo(result) {
    switch (result) {
        case 0:
            return {
                class: 'result-wild',
                text: 'Wild Type',
                interpretation: '정상형 (돌연변이 없음)'
            };
        case 1:
            return {
                class: 'result-heterozygous',
                text: 'Heterozygous',
                interpretation: '이형접합체 (정상형 + 돌연변이형)'
            };
        case 2:
            return {
                class: 'result-mutant',
                text: 'Mutant',
                interpretation: '돌연변이형 (동형접합체)'
            };
        default:
            return {
                class: 'result-no-call',
                text: 'No Call',
                interpretation: '결과 판정 불가'
            };
    }
}

// 시각화 결과 표시
function displayVisualizations(visualization1, visualization2) {
    // 패널 1 시각화
    if (visualization1) {
        visualizationPanel1.innerHTML = `
            <img src="data:image/png;base64,${visualization1}" 
                 class="img-fluid border rounded shadow" 
                 alt="패널 1 분석 결과"
                 style="max-width: 100%; height: auto;">
        `;
    } else {
        visualizationPanel1.innerHTML = `
            <div class="text-muted py-4">
                <i class="fas fa-exclamation-circle fa-2x mb-2"></i>
                <p>패널 1 시각화를 생성할 수 없습니다.</p>
            </div>
        `;
    }
    
    // 패널 2 시각화
    if (visualization2) {
        visualizationPanel2.innerHTML = `
            <img src="data:image/png;base64,${visualization2}" 
                 class="img-fluid border rounded shadow" 
                 alt="패널 2 분석 결과"
                 style="max-width: 100%; height: auto;">
        `;
    } else {
        visualizationPanel2.innerHTML = `
            <div class="text-muted py-4">
                <i class="fas fa-exclamation-circle fa-2x mb-2"></i>
                <p>패널 2 시각화를 생성할 수 없습니다.</p>
            </div>
        `;
    }
}

// 파일 드래그 앤 드롭 기능
function setupDragAndDrop() {
    const fileInputs = document.querySelectorAll('input[type="file"]');
    
    fileInputs.forEach(input => {
        const label = input.parentElement;
        
        // 드래그 오버 효과
        label.addEventListener('dragover', function(e) {
            e.preventDefault();
            label.classList.add('drag-over');
        });
        
        // 드래그 리브 효과
        label.addEventListener('dragleave', function(e) {
            e.preventDefault();
            label.classList.remove('drag-over');
        });
        
        // 드롭 처리
        label.addEventListener('drop', function(e) {
            e.preventDefault();
            label.classList.remove('drag-over');
            
            const files = e.dataTransfer.files;
            if (files.length > 0) {
                input.files = files;
                // 파일 선택 이벤트 트리거
                input.dispatchEvent(new Event('change', { bubbles: true }));
            }
        });
    });
}

// 파일 선택 시 파일명 표시
function setupFileDisplay() {
    const fileInputs = document.querySelectorAll('input[type="file"]');
    
    fileInputs.forEach(input => {
        input.addEventListener('change', function() {
            const fileName = this.files[0]?.name;
            const formText = this.nextElementSibling;
            
            if (fileName) {
                formText.innerHTML = `
                    <i class="fas fa-check-circle text-success me-1"></i>
                    선택된 파일: <strong>${fileName}</strong>
                `;
            } else {
                formText.innerHTML = 'ABI 또는 FSA 파일을 선택하세요.';
            }
        });
    });
}

// 페이지 로드 시 초기화
document.addEventListener('DOMContentLoaded', function() {
    setupDragAndDrop();
    setupFileDisplay();
    
    // 폼 리셋 기능
    uploadForm.addEventListener('reset', function() {
        hideError();
        hideResults();
        setLoadingState(false);
    });
});

// 키보드 단축키
document.addEventListener('keydown', function(e) {
    // Ctrl/Cmd + Enter로 분석 시작
    if ((e.ctrlKey || e.metaKey) && e.key === 'Enter') {
        if (!analyzeBtn.disabled) {
            uploadForm.dispatchEvent(new Event('submit'));
        }
    }
});

// 브라우저 뒤로가기/앞으로가기 처리
window.addEventListener('popstate', function(e) {
    if (e.state && e.state.results) {
        displayResults(e.state.results);
    }
});

// 결과를 브라우저 히스토리에 저장
function saveResultsToHistory(results) {
    history.pushState({ results: results }, '', window.location.href);
}
