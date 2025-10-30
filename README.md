# AutoCalling Web Application

유전자 분석을 위한 웹 기반 자동 호출 시스템입니다. ABI/FSA 파일을 업로드하여 다이플로타입 분석을 수행할 수 있습니다.

## 🚀 주요 기능

- **다중 패널 분석**: 2개의 컨트롤 패널과 2개의 테스트 패널을 동시에 분석
- **Ratio 기반 다중 분석**: 3가지 다른 ratio (0.50, 0.39, 0.28)로 분석하여 정확도 향상
- **실시간 시각화**: 분석 결과를 그래프로 시각화
- **다이플로타입 후보 제시**: 분석 결과에 따른 유전자형 후보 제시
- **결과 일관성 검증**: 여러 ratio에서의 결과 일관성 확인

## 📋 요구사항

- Python 3.7+
- Flask
- BioPython
- NumPy
- SciPy
- Matplotlib
- Pandas

## 🛠️ 설치 및 실행

### 1. 저장소 클론
```bash
git clone https://github.com/junfn1018-cmd/autocalling_web.git
cd autocalling_web
```

### 2. 가상환경 생성 및 활성화
```bash
python -m venv venv
source venv/bin/activate  # Windows: venv\Scripts\activate
```

### 3. 의존성 설치
```bash
pip install -r requirements.txt
```

### 4. 애플리케이션 실행
```bash
python app.py
```

### 5. 웹 브라우저에서 접속
```
http://localhost:5004
```

## 📁 프로젝트 구조

```
autocalling_web/
├── app.py                 # Flask 웹 애플리케이션
├── autocalling_web.py     # 핵심 분석 로직
├── templates/
│   └── index.html         # 웹 인터페이스
├── static/
│   ├── css/
│   │   └── style.css      # 스타일시트
│   └── js/
│       └── main.js        # JavaScript 로직
├── uploads/               # 업로드된 파일 저장소
├── requirements.txt       # Python 의존성
└── README.md             # 프로젝트 문서
```

## 🔬 사용 방법

1. **파일 업로드**: 4개의 FSA/ABI 파일을 업로드
   - 컨트롤 패널 1 (S1)
   - 컨트롤 패널 2 (S2)  
   - 테스트 패널 1 (S1)
   - 테스트 패널 2 (S2)

2. **분석 실행**: "분석 시작" 버튼 클릭

3. **결과 확인**: 
   - 다이플로타입 후보
   - 상세 분석 결과
   - 시각화 그래프

## 📊 분석 과정

1. **윈도우 생성**: 각 패널별 분석 윈도우 생성
2. **피크 탐지**: 컨트롤 샘플에서 기준 피크 탐지
3. **윈도우 수정**: 탐지된 피크를 기반으로 윈도우 조정
4. **다중 분석**: 3가지 ratio로 반복 분석
5. **결과 통합**: 패널별 결과를 결합하여 최종 다이플로타입 결정

## 🧬 지원 파일 형식

- **ABI 파일**: Applied Biosystems 형식
- **FSA 파일**: Applied Biosystems 형식

## ⚠️ 주의사항

- 업로드된 파일은 `uploads/` 폴더에 임시 저장됩니다
- 분석 완료 후 파일을 정리하는 것을 권장합니다
- 대용량 파일의 경우 분석 시간이 오래 걸릴 수 있습니다

## 🤝 기여하기

1. Fork the Project
2. Create your Feature Branch (`git checkout -b feature/AmazingFeature`)
3. Commit your Changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the Branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

## 📄 라이선스

이 프로젝트는 MIT 라이선스 하에 배포됩니다. 자세한 내용은 `LICENSE` 파일을 참조하세요.

## 📞 문의

프로젝트에 대한 문의사항이 있으시면 이슈를 생성해 주세요.

---

**개발자**: [Kim]  
**버전**: 1.0.0  
**최종 업데이트**: 2025년 10월
