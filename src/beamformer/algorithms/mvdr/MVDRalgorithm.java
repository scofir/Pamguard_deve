/*
 *  PAMGUARD - Passive Acoustic Monitoring GUARDianship.
 * To assist in the Detection Classification and Localisation
 * of marine mammals (cetaceans).
 *
 * Copyright (C) 2006
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 */



package beamformer.algorithms.mvdr;

import org.apache.commons.math3.complex.Complex; // NOTE: references to Complex are for this class, not the fftManager.Complex class
import org.apache.commons.math3.complex.ComplexField;
import org.apache.commons.math3.linear.ArrayFieldVector;
import org.apache.commons.math3.linear.FieldDecompositionSolver;
import org.apache.commons.math3.linear.FieldLUDecomposition;
import org.apache.commons.math3.linear.FieldMatrix;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.SingularMatrixException;

import Acquisition.AcquisitionProcess;
import Array.ArrayManager;
import Array.PamArray;
import Localiser.algorithms.PeakPosition;
import Localiser.algorithms.PeakSearch;
import PamController.PamController;
import PamDetection.LocContents;
import PamUtils.PamCalendar;
import PamUtils.PamUtils;
import PamUtils.complex.ComplexArray;
import PamView.dialog.warn.WarnOnce;
import beamformer.BeamAlgorithmParams;
import beamformer.BeamFormerBaseProcess;
import beamformer.algorithms.BeamFormerAlgorithm;
import beamformer.algorithms.BeamInformation;
import beamformer.continuous.BeamFormerDataBlock;
import beamformer.continuous.BeamFormerDataUnit;
import beamformer.continuous.BeamFormerProcess;
import beamformer.continuous.BeamOGramDataBlock;
import beamformer.continuous.BeamOGramDataUnit;
import beamformer.loc.BeamFormerLocalisation;
import beamformer.plot.BeamOGramLocalisation;
import fftManager.FFTDataUnit;
import pamMaths.PamVector;

/**
 * The algorithm class for the Minimum Variance Distortionless Response beamformer.  
 * 
 * @author mo55
 *
 */
public class MVDRalgorithm implements BeamFormerAlgorithm {

	/**
	 * The provider creating this algorithm
	 */
	private MVDRProvider mvdrProvider;
	
	/**
	 * Link back to the Beamformer process running the show
	 */
	private BeamFormerBaseProcess beamProcess;
	
	/**
	 * The parameters to use
	 */
	private MVDRParams mvdrParams;
	
	/**
	 * The sequence number to start with when creating beams
	 */
	private int firstBeamNum;
	
	/** 
	 * Output datablock, for all beams in all groups.
	 */
	private BeamFormerDataBlock beamformerOutput;

	/**
	 * The sequence number to use when creating a beamogram
	 */
	private int beamogramNum;

	/** 
	 * Output datablock, for all beams in all groups.
	 */
	private BeamOGramDataBlock beamOGramOutput;

	/**
	 * boolean indicating whether we are processing individual beams (true) or not (false)
	 */
	private boolean thisHasBeams=false;
	
	/**
	 * boolean indicating whether this is a BeamOGram (true) or not (false)
	 */
	private boolean thisHasABeamOGram=false;
	
	/**
	 * An array of PamVector objects describing the desired beam directions.  This is only used for individual beam
	 * analysis, not BeamOGrams.  It is initialised as empty, so that any calls to .length() will return 0 instead of a
	 * NullPointerException.
	 */
	private PamVector[] beamDirs = new PamVector[0];
	
	private BeamFormerLocalisation[] beamLocalisations;

	/**
	 * Object array containing beam information related to the beamogram (1 beam / look direction).  The first
	 * array is the secondary angle (aka slant, for a linear horizontal array), the angle relative to perpendicular
	 * to the primary array axis.  The second array is the main angle, in the direction of the primary array axis
	 */
	private MVDRBeam[][] beamogramBeams;

	/**
	 * A sequence map for the beamogram.  This is similar to a channelmap, but there is only one bit set.
	 */
	private int beamogramSeqMap;

	/**
	 * Object array containing individual beam information.
	 */
	private MVDRBeam[] individBeams;

	/**
	 * A sequence map for the beams.  This is similar to a channelmap, in that each bit is set for a certain beam.
	 */
	private int sequenceMap;

	/**
	 * The order of channels in the incoming FFT data units array
	 */
	protected int[] chanOrder = null;


	/**
	 * Keep frequency info in the beamOGram.
	 */
	private boolean keepFrequencyInfo = false;
	
	/**
	 * The inverse of the CSDM.  Each array index holds the CSDM for an fft bin.
	 */
	private FieldMatrix<Complex>[] Rinv;
	
	/**
	 * List of all channels in this group (derived from the channel map found in the parameters object)
	 */
	int[] channelList;

	/**
	 * The index in the fft data unit to start processing at.  Used only for the beamogram, and can be changed dynamically
	 * during processing
	 */
	private int beamOStartBin;

	/**
	 * The index in the fft data unit to end processing at.  Used only for the beamogram, and can be changed dynamically
	 * during processing
	 */
	private int beamOEndBin;

	/**
	 * Number of threads to process new data on
	 */
	int nBeamOThreads = 2;
	
	/**
	 * List of the beamformer threads
	 */
	private Thread[] beamThreads = new Thread[nBeamOThreads];
	
	private PeakSearch peakSearch = new PeakSearch(true);

	/**
	 * @param beamogramNum 
	 * @param firstSeqNum 
	 * @param mvdrParams 
	 * @param beamFormerProcess 
	 * @param mvdrProvider 
	 * 
	 */
	public MVDRalgorithm(MVDRProvider mvdrProvider, BeamFormerBaseProcess beamFormerProcess, MVDRParams mvdrParams, int firstSeqNum, int beamogramNum) {
		super();
		this.mvdrProvider = mvdrProvider;
		this.beamProcess = beamFormerProcess;
		this.mvdrParams = mvdrParams;
		this.firstBeamNum = firstSeqNum;
		this.beamogramNum = beamogramNum;
		this.beamformerOutput = beamFormerProcess.getBeamFormerOutput();
		this.beamOGramOutput = beamFormerProcess.getBeamOGramOutput();
		
		// create the vector of beam directions for individual beams
		if (mvdrParams.getNumBeams()>0) {
			int locContents = LocContents.HAS_BEARING | LocContents.HAS_AMBIGUITY | LocContents.HAS_BEARINGERROR;
			thisHasBeams = true;
			beamDirs = new PamVector[mvdrParams.getNumBeams()]; // 创建存储单独波束方向的 PamVector 数组
			int[] headings = mvdrParams.getHeadings(); // 获取波束的方向角数组
			int[] slants = mvdrParams.getSlants(); // 获取波束的倾斜角数组
			beamLocalisations = new BeamFormerLocalisation[mvdrParams.getNumBeams()]; // 创建存储单独波束位置信息的 BeamFormerLocalisation 数组
			for (int i=0; i<mvdrParams.getNumBeams(); i++) {
				beamDirs[i]=PamVector.fromHeadAndSlant(headings[i], slants[i]); // 使用方向角和倾斜角创建 PamVector 对象
				double beamErr = 180.;
				if (mvdrParams.getNumBeams() > 1) {
					if (i == 0 || i == mvdrParams.getNumBeams()-1) {
						beamErr = Math.abs(headings[1]-headings[0]);
					}
					else {
						beamErr = Math.abs(headings[i+1]-headings[i-1]) / 2.; // 计算中间波束的方向角之差的一半作为波束误差
					}
				}
//				beamLocalisations[i] = new BeamFormerLocalisation(null, 
//						locContents, 
//						mvdrParams.getChannelMap(), 
//						Math.toRadians(headings[i]), 
//						Math.toRadians(beamErr));
			}
		}
		
		// get the parameters for the beamogram
		if (mvdrParams.getNumBeamogram()>0) {
			thisHasABeamOGram = true;
		}
	}

	/* (non-Javadoc)
	 * @see beamformer.algorithms.BeamFormerAlgorithm#prepare()
	 */
	@Override
	public void prepare() {	// Prepare the algorithm. Gets called just before it starts. It's here that steering vectors, etc. should get calculated. 
		// 计算导向矢量
		// get a list of all channels in this channel map, and the corresponding element locations
		channelList = PamUtils.getChannelArray(mvdrParams.getChannelMap()); // 获取通道映射中的所有通道列表和对应的元素位置
		AcquisitionProcess sourceProcess = null;
		try {
//			sourceProcess = (AcquisitionProcess) beamProcess.getSourceProcess();
			sourceProcess = (AcquisitionProcess) beamProcess.getFftDataSource().getSourceProcess(); // 尝试获取采集进程
		}
		catch (ClassCastException e) {
			String title = "Error finding Acquisition module";
			String msg = "There was an error trying to find the Acquisition module.  " +
					"The beamformer needs this information in order to run.  There will be no output until " +
					"a valid Acquisition module is added and the Pamguard run is restarted.";
			String help = null;
			int ans = WarnOnce.showWarning(PamController.getMainFrame(), title, msg, WarnOnce.WARNING_MESSAGE, help, e);
			sourceProcess=null;
			e.printStackTrace();
			return;
		}
		if (sourceProcess==null) {
			String title = "Error finding Acquisition module";
			String msg = "There was an error trying to find the Acquisition module.  " +
					"The beamformer needs this information in order to run.  There will be no output until " +
					"a valid Acquisition module is added and the Pamguard run is restarted.";
			String help = null;
			int ans = WarnOnce.showWarning(PamController.getMainFrame(), title, msg, WarnOnce.WARNING_MESSAGE, help, null);
			return;
		}
		
		ArrayManager arrayManager = ArrayManager.getArrayManager();
		PamArray currentArray = arrayManager.getCurrentArray();
		PamVector[] elementLocs = new PamVector[channelList.length];
		long now = PamCalendar.getTimeInMillis();
		int hydrophoneMap = 0;
		for (int i=0; i<channelList.length; i++) {
			int hydrophone = sourceProcess.getAcquisitionControl().getChannelHydrophone(channelList[i]); // 获取通道对应的水听器编号
			elementLocs[i] = currentArray.getAbsHydrophoneVector(hydrophone, now); // 获取当前时刻下水听器的绝对位置
			hydrophoneMap |= 1<<hydrophone;         // 更新水听器映射
		}
		/*
		 * Do some normalisation of those vectors. Start by making everything relative
		 * to the average position. 
		 */
		PamVector arrayCenter = PamVector.mean(elementLocs);	// 计算元素位置的平均值，得到阵列中心点
		for (int i = 0; i < channelList.length; i++) {
			elementLocs[i] = elementLocs[i].sub(arrayCenter);	// 将每个元素位置向量减去阵列中心点，使其相对于平均位置
		}
		/*
		 * Now get the principle axis vector of the array. 
		 */
		int arrayShape = arrayManager.getArrayType(hydrophoneMap); // 2 for  aline array. 
		PamVector[] arrayAxis = arrayManager.getArrayDirection(hydrophoneMap); // will be single vector for a line array. 获取阵列的主轴向量，对于线阵列来说，返回值是一个单一的向量
		
		// Create the beams for the beamogram.  Step through the look angles (based on the beamOGramAngles variable) and create a vector for each.  Use
		// the same sequence number, since when processing all beams will be added together anyway.  For the element weights and frequency range, use
		// the last index in the weights and freqRange variables (the last index = number of individual beams)
		int[] boAngles = mvdrParams.getBeamOGramAngles(); // 获取波束图的观察角度数组 {0，180，2步长}
		if (thisHasABeamOGram && boAngles != null && boAngles.length >= 2) {
			int nMainBeams = (int) ((mvdrParams.getBeamOGramAngles()[1]-mvdrParams.getBeamOGramAngles()[0])/mvdrParams.getBeamOGramAngles()[2]+1); // 主方向波束数量 （180-0）/2 +1=91
			int nSecBeams = (int) ((mvdrParams.getBeamOGramSlants()[1]-mvdrParams.getBeamOGramSlants()[0])/mvdrParams.getBeamOGramSlants()[2]+1); // 次方向（斜角）波束数量 0
			beamogramBeams = new MVDRBeam[nSecBeams][nMainBeams]; // 创建一个大小为[nSecBeams][nMainBeams]的MVDRBeam数组，用于存储波束
			beamogramSeqMap = PamUtils.SetBit(0, beamogramNum, true);
//			int arrayIdx = 0;
//			if (beamDirs != null) {
//				arrayIdx = beamDirs.length;
//			}
			PamVector beamDir = new PamVector();
			for (int j=0; j<nSecBeams; j++) {
				for (int i=0; i<nMainBeams; i++) {
					beamDir=PamVector.fromHeadAndSlant(mvdrParams.getBeamOGramAngles()[0]+i*mvdrParams.getBeamOGramAngles()[2],
							mvdrParams.getBeamOGramSlants()[0]+j*mvdrParams.getBeamOGramSlants()[2]);
					beamogramBeams[j][i] = new MVDRBeam(this, 
							mvdrParams.getChannelMap(), 
							beamogramNum, 
							beamDir, 
							elementLocs, 
							mvdrParams.getBeamOGramFreqRange(), 
							currentArray.getSpeedOfSound());
				} // 计算各个角度的导向矢量 beamogramBeams
			}
			double[] freqBins = mvdrParams.getBeamOGramFreqRange();     // 获取波束图的频率范围
			beamOStartBin = beamProcess.frequencyToBin(freqBins[0]);
			beamOEndBin = beamProcess.frequencyToBin(freqBins[1]);
			
		// Create the individual beams.  Loop through the look vectors created in the constructor.  Give each beam a unique seequence number
		} // 创建单个波束。循环遍历在构造函数中创建的look向量，并为每个波束分配唯一的序列号。
		if (thisHasBeams) {
			individBeams = new MVDRBeam[beamDirs.length];
			for (int i=0; i<beamDirs.length; i++) {
				sequenceMap = PamUtils.SetBit(0, firstBeamNum+i, true);
				individBeams[i] = new MVDRBeam(this, 
						mvdrParams.getChannelMap(), 
						firstBeamNum+i, 
						beamDirs[i], 
						elementLocs, 
						mvdrParams.getFreqRange()[i], 
						currentArray.getSpeedOfSound());
			}
		}
		
		// reset the channel order array
		this.clearChannelOrderList();
	}

	/**
	 * Process an array of FFTDataUnits.  The size of the array will equal the number of channels in this beam group
	 */
	@Override
	public void process(FFTDataUnit[] fftDataUnits) {
		
		// if this is a BeamOGram, we need to loop through all of the beams, averaging the magnitudes of the returned complex numbers and then
		// saving that to a new index in a ComplexArray.  This ComplexArray will then be used to create a new BeamFormerDataUnit.  Note that we
		// store the values into the ComplexArray in reverse order (from last bin to first), so that on the spectrogram display the value related
		// to heading=0deg (straight ahead) is at the top of the display, and the value related to heading=180deg (behind the boat) is at the
		// bottom.
		if (thisHasABeamOGram && beamogramBeams != null) {
			
			// calculate the inverse CSDM for the data over the currently-specified frequency bin range
//			this.prepNewData(fftDataUnits, beamogramBeams[0][0].getStartIdx(), beamogramBeams[0][0].getNumFFTBins()); this one calculates CSDM for entire range
			this.prepNewData(fftDataUnits, beamOStartBin, beamOEndBin-beamOStartBin); // 求解每个频点，信号协方差矩阵的逆 Rinv 维度（fftBin个ch*ch矩阵）
			
			// loop through the angles one beam at a time, processing the data and averaging the results
			int numAnlgeBins = beamogramBeams[0].length;//beamProcess.getFftDataSource().getFftLength()/2; 获取主角度数
			int numSlantBins = beamogramBeams.length; // 获取次角度数
			int nAllFBins = fftDataUnits[0].getFftData().length(); // 获取总频点数
			int nFBinsToKeep = keepFrequencyInfo ? nAllFBins : 1;
			double[][][] beamData = new double[nFBinsToKeep][numSlantBins][numAnlgeBins];
			/*
			 * Keep the order simple here - for the data, not the displays ! So bin 0 in the output 
			 * is bin 0 in angle. Don't reverse it.  DG. 27.07.17
			 * 
			 * Multithreading. Make multiple threads for this next part, looping 
			 * over the second (final) dimension, which is generaly the one having the  多线程处理。为这部分代码创建多个线程，循环遍历第二个（最后一个）维度，通常这是具有最多独立波束的维度。
			 * most separate beams. 
			 */
			if (nBeamOThreads > 1) {
				for (int i = 0; i < nBeamOThreads; i++) {
					beamThreads[i] = new Thread(new BeamOThread(fftDataUnits, beamData, i, nBeamOThreads, nAllFBins)); // 多线程计算每个方位角（二维则包括方位和俯仰角） 功率值
					beamThreads[i].start();
				}
				try {
					for (int i = 0; i < nBeamOThreads; i++) {
						beamThreads[i].join();
					}
				} catch (InterruptedException e) {
					e.printStackTrace();
				}

				
			}
			else {
				for (int j=0; j<beamogramBeams.length; j++) {
					for (int i=0; i<beamogramBeams[0].length; i++) {
						//					ComplexArray summedData = beamogramBeams[j][i].process(Rinv);
						ComplexArray summedData = beamogramBeams[j][i].process(Rinv, beamOStartBin, beamOEndBin-beamOStartBin,nAllFBins); // 计算每个方位角（二维则包括方位和俯仰角）全频点的功率谱（复数向量）
						/*
						 * Do a search for NaN values here ...
						 */
						double[] d = summedData.getData();
//						int nNan = 0;
//						for (int r = 1; r < d.length; r++) {
//							if (Double.isNaN(d[r])) {
//								nNan ++;
//							}
//						}
//						if (nNan > 0) {
//							summedData = beamogramBeams[j][i].process(Rinv, beamOStartBin, beamOEndBin-beamOStartBin,nAllFBins);
//						}
						//				double summedMag = 0.;
						//				for (int k=0; k<summedData.length(); k++) {
						//					summedMag += summedData.magsq(k);
						//				}
						//					for (int k=0; k<summedData.length(); k++) {
						for (int k=beamOStartBin; k<beamOEndBin; k++) {
							if (keepFrequencyInfo) {
								beamData[k][j][i] = summedData.magsq(k);
							}
							else {
								beamData[0][j][i] += summedData.magsq(k); //计算每个方位角（二维则包括方位和俯仰角）的功率值（实数）
							}
						}
					}
				}
			}
			
			// create a new data unit
			BeamOGramDataUnit newUnit = new BeamOGramDataUnit(fftDataUnits[0].getTimeMilliseconds(),
					mvdrParams.getChannelMap(),
					beamogramSeqMap, 
					fftDataUnits[0].getStartSample(),
					fftDataUnits[0].getSampleDuration(),
					beamData, 
					fftDataUnits[0].getFftSlice());
			
			// now get the best angle and add it as a localisation. 
			double meanLev = 0;
			double[] aveData = newUnit.getAngle1Data(true); // 输出avedata,是一维的方位角功率值，bins和俯仰角维度都被平均了
			for (int i = 0; i < aveData.length; i++) {
				meanLev += aveData[i];
			}
			meanLev /= aveData.length;
			PeakPosition peakPos = peakSearch.interpolatedPeakSearch(aveData); // 使用峰值搜索算法查找峰值位置 针对一维向量（方位角的功率值）
			if (peakPos.getHeight()/meanLev > 2) {                             // 如果峰值高度相对于平均值超过2倍，则认为是有效的峰值
				int[] boga = mvdrParams.getBeamOGramAngles(); // 获取BeamOGram 方位角 {0， 180， 2（步长）}
				double bestAng = peakPos.getBin() * (double) boga[2] + (double) boga[0]; // 计算最佳方位角度
				newUnit.setLocalisation(new BeamOGramLocalisation(newUnit, mvdrParams.getChannelMap(), bestAng));
			}
			
			beamOGramOutput.addPamData(newUnit);
		}
		
		// if these are individual beams, loop through them one at a time and create a new BeamFormerDataUnit for each 如果这些是单独的波束，逐个循环它们并为每个创建一个新的BeamFormerDataUnit
		if (thisHasBeams) {
			int iBeam = 0;
			
			// set the counters and calculate the inverse CSDM for the first beam 为第一个波束计算信号协方差矩阵
			int lastStartIdx = individBeams[0].getStartIdx();
			int lastNumFFTBins = individBeams[0].getNumFFTBins();
			this.prepNewData(fftDataUnits, individBeams[0].getStartIdx(), individBeams[0].getNumFFTBins());
			
			// loop through the beams one at a time, saving each result to the beamformer output data block 设置计数器并计算第一个波束的逆CSDM（交叉谱密度矩阵）
			for (MVDRBeam beam : individBeams) {
				
				// if this beam has a different starting index or number of fft bins as the previous, recalculate the inverse CSDM
				if (beam.getStartIdx()!=lastStartIdx || beam.getNumFFTBins()!=lastNumFFTBins) {
					this.prepNewData(fftDataUnits, beam.getStartIdx(), beam.getNumFFTBins());
					lastStartIdx = beam.getStartIdx();
					lastNumFFTBins = beam.getNumFFTBins();
				}
				ComplexArray summedData = beam.process(Rinv);
				BeamFormerDataUnit newUnit = new BeamFormerDataUnit(fftDataUnits[0].getTimeMilliseconds(),
						mvdrParams.getChannelMap(),
						beam.getSequenceMap(), 
						fftDataUnits[0].getStartSample(),
						fftDataUnits[0].getSampleDuration(),
						summedData, 
						fftDataUnits[0].getFftSlice());
//				System.out.println("Beam " + String.valueOf(beam.getSequenceMap()) + " at 150Hz = " +String.valueOf(summedData.getReal(19)) + "+" + String.valueOf(summedData.getImag(19))+"i");
				newUnit.setLocalisation(beamLocalisations[iBeam]);
				beamformerOutput.addPamData(newUnit);
				iBeam++;
			}
		}
	}

	/**
	 * Multithreading inner class to process new data
	 * @author mo55
	 *
	 */
	private class BeamOThread implements Runnable {

		private FFTDataUnit[] fftDataUnits;
		private double[][][] beamData;
		private int firstBeamIndex;
		private int beamIndexStep;
		private int nFBins;

		public BeamOThread(FFTDataUnit[] fftDataUnits, double[][][] beamData, int firstBeamIndex, int beamIndexStep, int nFBins) {
			super();
			this.fftDataUnits = fftDataUnits;
			this.beamData = beamData;
			this.firstBeamIndex = firstBeamIndex;
			this.beamIndexStep = beamIndexStep;
			this.nFBins = nFBins;
		}

		@Override
		public void run() {
			for (int j=0; j<beamogramBeams.length; j++) {
				for (int i=firstBeamIndex; i<beamogramBeams[0].length; i+=beamIndexStep) {
					ComplexArray summedData = beamogramBeams[j][i].process(Rinv, beamOStartBin, beamOEndBin-beamOStartBin,nFBins); // 计算MVDR的功率谱
					for (int k=beamOStartBin; k<beamOEndBin; k++) {
						if (keepFrequencyInfo) {
							beamData[k][j][i] = summedData.magsq(k);
						}
						else {
							beamData[0][j][i] += summedData.magsq(k); // 计算功率值（每个角度每个频点单独）的模
						}
					}
					//				aveData.set(j, summedMag/summedData.length(), 0.);
					//				aveData[j] = summedMag/summedData.length();
				}
			}
			
		}
		
	}

	/**
	 * Calculate the inverse CSDM (Rinv field) based on the current fftDataUnits. 循环每个频点，求该频点处信号协方差矩阵的逆 Rinv （fftBin*ch*ch）
	 * 
	 * @param fftDataUnits the new data to process
	 * @param startIdx the fft bin to start the processing at
	 * @param numFFTBins the number of FFT bins to process
	 */
	private void prepNewData(FFTDataUnit[] fftDataUnits, int startIdx, int numFFTBins) {
		Rinv = (FieldMatrix<Complex>[]) new FieldMatrix[numFFTBins]; // cannot create an array with generics, so need to do this awkward cast instead
		
		// if we don't know the order of channels for the incoming FFT data units, determine that now
		if (chanOrder==null) {
			getChannelOrder(fftDataUnits);
		}

		// loop over the FFT bins calculating an inverse CSDM for each one
		for (int fftBin=0; fftBin<numFFTBins; fftBin++) {
				
			// compile the FFT values for each hydrophone in this FFT bin into both a regular vector as well as a
			// conjugate vector, so that later we can calculate the CSDM (outer product) easily.  Make sure to compile them
			// in the order that matches the steering vector order.
			// Note: the Complex objects used here are from the Apache Commons Math library org.apache.commons.math3.complex.Complex,
			// not the fftManager.Complex class
			ArrayFieldVector<Complex> dataVec = new ArrayFieldVector<Complex>(ComplexField.getInstance(), fftDataUnits.length); // fftDataUnits.length 可以理解为通道数
			ArrayFieldVector<Complex> dataVecConj = new ArrayFieldVector<Complex>(ComplexField.getInstance(), fftDataUnits.length);
			for (int chan=0; chan<fftDataUnits.length; chan++) {
				dataVec.setEntry(chan, new Complex(fftDataUnits[chanOrder[chan]].getFftData().get(startIdx+fftBin).real, fftDataUnits[chanOrder[chan]].getFftData().get(startIdx+fftBin).imag));
				dataVecConj.setEntry(chan, new Complex(fftDataUnits[chanOrder[chan]].getFftData().get(startIdx+fftBin).real, -1*fftDataUnits[chanOrder[chan]].getFftData().get(startIdx+fftBin).imag));
			}
			
			// calculate the CSDM of the FFT data
			FieldMatrix<Complex> R = dataVec.outerProduct(dataVecConj); // 计算的CSDM 是协方差矩阵（Cross-Spectral Density Matrix），outerProduct计算外积 （1*ch） x （1*ch） = （ch*ch）=R
			
			// add diagonal loading to the matrix.  First calculate the trace value (the sum of the diagonal components)
			// and divide that by the square of the number of channels.  Add that to the diagonal elements in the matrix R
			Complex traceVal = R.getTrace(); // 对协方差矩阵，对角加载操作，对角元素的加载量为 迹值（对角线元素的总和）除以通道数的平方。
			Complex noiseVal = traceVal.divide(Math.pow(fftDataUnits.length,2));
			for (int chan=0; chan<fftDataUnits.length; chan++) {
				R.addToEntry(chan, chan, noiseVal);
			}
			
			// divide by the number of channels.  This was done in the Matlab program.  Not really necessary, but it helps to bring
			// the data into a scale closer to the original FFT data.
//			R.scalarMultiply(new Complex(1./fftDataUnits.length));
			
			
			// calculate the inverse matrix for this fft bin
			FieldDecompositionSolver<Complex> solver = new FieldLUDecomposition<Complex>(R).getSolver(); // 通过getSolver() 方法获取到 LU 分解对象的求解器，用于后续矩阵求逆
			FieldMatrix<Complex> singleRinv = MatrixUtils.createFieldMatrix(ComplexField.getInstance(), R.getRowDimension(), R.getColumnDimension()); // 创建一个空的复数矩阵存储逆矩阵运算结果 
			try {
				singleRinv = solver.getInverse();

			} catch (SingularMatrixException ex) {
			}
			
			// save the inverse CSDM to the field array
			Rinv[fftBin] = singleRinv;
		}
	}

	/**
	 * loop over the number of channels and determine the order of the hydrophones in the FFTDataUnits object.
	 * Create a look up table to match the order of hydrophones in the fftDataUnits array to the order
	 * in the steeringVecs array
	 * 
	 * @param fftDataUnits
	 */
	protected void getChannelOrder(FFTDataUnit[] fftDataUnits) {
		chanOrder = new int[channelList.length];
		for (int i = 0; i < fftDataUnits.length; i++) {
			int chanToMatch = PamUtils.getSingleChannel(fftDataUnits[i].getChannelBitmap());
			for (int j=0; j<channelList.length; j++) {
				if (channelList[j]==chanToMatch) {
					chanOrder[j]=i;
					break;
				}
			}
		}
	}

	/**
	 * clear the list of channel orders
	 */
	public void clearChannelOrderList() {
		chanOrder=null;
	}

	/* (non-Javadoc)
	 * @see beamformer.algorithms.BeamFormerAlgorithm#getNumBeams()
	 */
	@Override
	public int getNumBeams() {
		return beamDirs.length;
	}

	/* (non-Javadoc)
	 * @see beamformer.algorithms.BeamFormerAlgorithm#getNumBeamogramAngles()
	 */
	@Override
	public int getNumBeamogramAngles() {
		int[] boAngles = mvdrParams.getBeamOGramAngles();
		if (boAngles == null) return 1;
		int numAngles = (int) ((mvdrParams.getBeamOGramAngles()[1]-mvdrParams.getBeamOGramAngles()[0])/mvdrParams.getBeamOGramAngles()[2]+1);
		return numAngles;
	}

	/* (non-Javadoc)
	 * @see beamformer.algorithms.BeamFormerAlgorithm#getBeamInformation(int)
	 */
	@Override
	public BeamInformation getBeamInformation(int iBeam) {
		return null;
	}

	/* (non-Javadoc)
	 * @see beamformer.algorithms.BeamFormerAlgorithm#thereIsABeamogram()
	 */
	@Override
	public boolean thereIsABeamogram() {
		return (thisHasABeamOGram);
	}
	
	/**
	 * @return the beamProcess
	 */
	public BeamFormerBaseProcess getBeamProcess() {
		return beamProcess;
	}

	@Override
	public void setKeepFrequencyInformation(boolean keep) {
		keepFrequencyInfo  = true;
	}

	@Override
	public void setFrequencyBinRange(int binFrom, int binTo) {
		beamOStartBin = binFrom;
		beamOEndBin = binTo;		
	}

}
