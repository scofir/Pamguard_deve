package Localiser.algorithms.timeDelayLocalisers.bearingLoc;

import Array.ArrayManager;
import Array.PamArray;
import Jama.LUDecomposition;
import Jama.Matrix;
import Jama.QRDecomposition;
import PamUtils.PamUtils;
import pamMaths.PamVector;

public class LSQBearingLocaliser implements BearingLocaliser {

	private int hydrophoneBitMap;
	private long timeMillis;
	private double timingError;
	
	private Matrix weightedHydrophoneVectors;          // 加权水听器向量矩阵
	private Matrix hydrophoneVectors;                  // 水听器向量矩阵
	private Matrix hydrophoneErrorVectors;             // 水听器误差向量矩阵
	private double[] hydrophoneSpacing;                // 水听器间距数组
	private PamArray currentArray;                     // 当前阵列
	private int arrayType;                             // 阵列类型
	private PamVector[] arrayAxis;                     // 阵列轴
	private QRDecomposition qrHydrophones;             // 水听器的QR分解对象
	private double[] fitWeights;                       // 拟合权重数组

	public LSQBearingLocaliser(int hydrophoneBitMap, long timeMillis, double timingError) {
		this.hydrophoneBitMap = hydrophoneBitMap;
		this.timeMillis = timeMillis;
		this.timingError = timingError;
	}

	@Override
	public void prepare(int[] arrayElements, long timeMillis, double timingError) {
		/*
		 * Set up the matrixes of inter hydrophone vectors. 计算 加权水听器向量矩阵，水听器向量矩阵，水听器误差向量矩阵 维度都是（Npairs，3）用于表示，每对水听器（任意2个）间距的信息
		 */
		this.timingError = timingError;

		hydrophoneBitMap = PamUtils.makeChannelMap(arrayElements); // 生成水听器位图
		ArrayManager arrayManager = ArrayManager.getArrayManager();
		currentArray = arrayManager.getCurrentArray();
		arrayType = arrayManager.getArrayShape(currentArray, hydrophoneBitMap);

		arrayAxis = arrayManager.getArrayDirections(currentArray, hydrophoneBitMap);
		
		int nHyd = arrayElements.length; // 水听器数量
		int nDelay = (nHyd*(nHyd-1))/2;  // 延迟值的数量
		weightedHydrophoneVectors = new Matrix(nDelay, 3); // 初始化矩阵和数组
		hydrophoneVectors = new Matrix(nDelay, 3);
		hydrophoneErrorVectors = new Matrix(nDelay, 3);
		hydrophoneSpacing = new double[nDelay];
		fitWeights = new double[nDelay];
		double c = currentArray.getSpeedOfSound(); // 声速
		int iRow = 0; 							   // 行索引
		for (int i = 0; i < nHyd; i++) {
			PamVector vi = currentArray.getAbsHydrophoneVector(i, timeMillis); // Get the true location vector for  a hydrophone element  获取绝对位置的向量坐标
			for (int j = i+1; j <nHyd; j++) {
				PamVector vj = currentArray.getAbsHydrophoneVector(j, timeMillis);
				PamVector v = vj.sub(vi);                                             // 计算向量差
				hydrophoneSpacing[iRow] = v.norm();                                   // 计算向量差L2范数表征水听器间距 即两坐标的欧几里得距离
				PamVector uv = v.getUnitVector();									// 获取当前向量的单位向量 （当前向量的每个分量除以向量的长度）
				PamVector errorVec = currentArray.getSeparationErrorVector(i, j, timeMillis);    // 获取水听器i和j之间的误差向量
				/*
				 * Get the error component along the line of the pair and then calculate
				 * the weight as 1/the variance of the hydrophone separation. 
				 */
				double errorComponent = uv.dotProd(errorVec); // 误差向量与单位向量点积求和，得到误差向量在uv上的投影 可以用来量化误差向量在单位向量uv方向上的分量大小，正负号表示误差的方向，绝对值表示误差的大小。
				fitWeights[iRow] = Math.pow(v.norm()/errorComponent, 2); // 用于表示水听器之间距离的拟合权重 可以用来衡量误差向量errorVec在单位向量uv方向上对v范数的相对影响
				for (int e = 0; e < 3; e++) {
					weightedHydrophoneVectors.set(iRow, e, v.getElement(e)/c*fitWeights[iRow]); // 加权水听器向量矩阵的维度 （nDelay，3），计算并填充到矩阵，可理解为每个水听器对的距离延迟（无法方向信息）
					hydrophoneVectors.set(iRow, e, v.getElement(e)/c);
					hydrophoneErrorVectors.set(iRow, e, errorVec.getElement(e)/c);
//					hydrophoneUnitVectors.set(iRow, e, uv.getElement(e));
				}
				iRow++;
			}
		}
//		luHydrophoneUnitMatrix = new LUDecomposition(hydrophoneUnitVectors);
		qrHydrophones = new QRDecomposition(weightedHydrophoneVectors); // 对加权水听器向量矩阵的QR分解对象
	}

	@Override
	public int getArrayType() {
		return arrayType;
	}

	@Override
	public int getHydrophoneMap() {
		return hydrophoneBitMap;
	}

	@Override
	public PamVector[] getArrayAxis() {
		return arrayAxis;
	}	 
	
	/* 
	 * @return true if a new grid needs to be created // 阵列位置改变时，需要更新prepare中的内容
	 */
	private boolean resetArray(long timeMillis){
		// 检查当前阵列是否为空或者时间发生了变化并且水听器定位器是可变的
		if (currentArray == null || (this.timeMillis!=timeMillis && currentArray.getHydrophoneLocator().isChangeable())){
			prepare(PamUtils.getChannelArray(hydrophoneBitMap), timeMillis, 1e-6);
			this.timeMillis = timeMillis;
			return true;
		}
			
		return false; 
	}

	@Override
	public double[][] localise(double[] delays, long timeMillis) { // 通过计算构建方程的最小二乘解，算出方位角和俯仰角
		resetArray(timeMillis);                                   // 重置阵列 更新阵列位置
//		qrHydrophones = new QRDecomposition(hydrophoneVectors);
		Matrix normDelays = new Matrix(delays.length, 1);
		for (int i = 0; i < delays.length; i++) {
			normDelays.set(i, 0, -delays[i]*fitWeights[i]); // 负号表示对延迟值进行反向修正。
		}
//		Matrix soln = luHydrophoneUnitMatrix.solve(normDelays);
		Matrix soln2 = qrHydrophones.solve(normDelays); // 解决线性方程组 A * X = B，其中 A 是进行了 QR 分解的矩阵qrHydrophones，B 是右侧向量 normDelays。 方程解X即soln2维度是（3，1）表示最小二乘拟合的方向坐标
		double[][] angs = new double[2][2];
		PamVector v = new PamVector(soln2.get(0, 0), soln2.get(1,0), soln2.get(2, 0)); // v是计算出的声源坐标，可归一化后用atan和asin求方位角和俯仰角
//		System.out.printf("Vector Norm = %4.3f: ", v.norm());
		double m = v.normalise();                // 归一化向量,消除长度的影响 此时的v已经被归一化，m只是用来存函数返回值，无意义
		angs[0][0] = Math.PI/2. - Math.atan2(v.getElement(0),v.getElement(1));     // 计算方位角，公式是从三维坐标计算角度
		angs[0][1] = Math.asin(v.getElement(2));                                   // 计算俯仰角
		
//		timingError = 1e-5;
		// now take a look at angle errors
		double oneDeg = Math.PI/180.;
		
		double testDeg = 5;
		double[][] er = new double[2][20];
		for (int i = 0; i < 20; i++) {
			testDeg = 1+i;
		// pick points testDeg degrees either side of angs and change one at a time 选择点测试角度两边的度数，每次改变1度
		double aDiff = testDeg * oneDeg;
		double a;
		double l1, l2, l3, l1a, l3a;
		l1 = logLikelihood(delays, angs[0][0] - aDiff, angs[0][1]); // 正负20度按步长1度，求解对数似然度
		l2 = logLikelihood(delays, angs[0][0], angs[0][1]);
		l3 = logLikelihood(delays, angs[0][0] + aDiff, angs[0][1]);
		er[0][i] = angs[1][0] = Math.sqrt(1./(l1+l3-2*l2))*aDiff;   // 计算方向角误差 （计算曲率）
		l1a = logLikelihood(delays, angs[0][0], angs[0][1] - aDiff);
//		l2 = logLikelihood(delays, angs[0][0], angs[0][1]);
		l3a = logLikelihood(delays, angs[0][0], angs[0][1] + aDiff);
		er[1][i] = angs[1][1] = Math.sqrt(1./(l1a+l3a-2*l2))*aDiff; // 计算俯仰角误差
		}
		
		
		
		
//		double ll[] = new double[21];
//		double a[] = new double[2];
//		timingError = 1.e-4;
//		a[1] = angs[0][1];
//		for (int i = 0; i < ll.length; i++) {
//			a[0] = angs[0][0] + (-10 + i)*oneDeg;
//			ll[i] = logLikelihood(delays, a);
//		}
		return angs; // angs [ 方位角 ，俯仰角
	}				 //		  方位角误差， 俯仰角误差]
	/**
	 * Log likelihood function
	 * @param delays time delays
	 * @param angle0 horizontal angle
	 * @param angle1 vertical angle
	 * @return log likelihood based on an estimate of errors. 
	 */
	public double logLikelihood(double[] delays, double angle0, double angle1) {
		Matrix whaleVec = new Matrix(3,1);
		whaleVec.set(0, 0, Math.cos(angle1)*Math.cos(angle0));
		whaleVec.set(1, 0, Math.cos(angle1)*Math.sin(angle0));
		whaleVec.set(2, 0, Math.sin(angle1));
		return logLikelihood(delays, whaleVec);
	}
	
	/**
	 * Calculate a log likelihood for a given pair of angles
	 * @param delays actual delays
	 * @param angles angles, horizontal and vertical. 
	 * @return
	 */
	public double logLikelihood(double[] delays, double[] angles) {
		return logLikelihood(delays, angles[0], angles[1]);
//		Matrix whaleVec = new Matrix(3,1);
//		whaleVec.set(0, 0, Math.cos(angles[1])*Math.cos(angles[0]));
//		whaleVec.set(1, 0, Math.cos(angles[1])*Math.sin(angles[0]));
//		whaleVec.set(2, 0, Math.sin(angles[1]));
//		return logLikelihood(delays, whaleVec);
	}
	/**
	 * Calculate a log likelihood for a given whale vector.  对数似然度函数常用于度量观测值与模型预测值之间的差异。
	 * @param delays actual delays
	 * @param whaleVector estimated whale vector. 
	 * @return
	 */
	public double logLikelihood(double[] delays, Matrix whaleVector) {
		Matrix times = hydrophoneVectors.times(whaleVector); // expected times for this whale position (note vecs are already divided by c) 计算鲸鱼位置对应的预期时间（注意向量已经除以c）
		Matrix timeErrors = hydrophoneErrorVectors.times(whaleVector);
		double c = currentArray.getSpeedOfSound();
		double dc = currentArray.getSpeedOfSoundError();
		Matrix timeErrors2 = times.times(dc/c/c);
		double chi = 0;
		for (int i = 0; i < times.getRowDimension(); i++) {
			double expectedVariance = Math.pow(timeErrors.get(i, 0),2) + Math.pow(timeErrors2.get(i, 0),2) + Math.pow(timingError, 2);  // 计算预期方差，包括时间误差平方、时间误差2平方和定时误差平方的总和
			chi = Math.pow((times.get(i, 0)+delays[i]), 2)/expectedVariance; // 计算卡方值，将预期时间与实际延迟时间之和的平方除以预期方差。chi 是通过将观测值与模型预测值之间的差异除以 expectedVariance 得到的标准化残差。
		}																	// 标准化残差是一种用于评估观测值与模型之间差异的指标，它表示观测值与模型预测值之间的差异相对于预期方差的大小。
		return chi/2;
	}

}
