package beamformer;

import Array.ArrayManager;
import Array.PamArray;
import PamUtils.PamUtils;
import beamformer.algorithms.BeamAlgorithmProvider;
import beamformer.algorithms.BeamFormerAlgorithm;
import fftManager.FFTDataUnit;
import pamMaths.PamVector;

/**
 * Handle channel groupings so that the actual algorithm
 * gets called with a complete array of all channels in 
 * it's group in one call. Do it here since every single 
 * algorithm will need to do this, so save them the effort. 
 * @author dg50
 *
 */
public class BeamGroupProcess {
	
	private BeamFormerBaseProcess beamFormerBaseProcess;
	private BeamAlgorithmProvider provider;
	private BeamFormerAlgorithm beamFormerAlgorithm;
	private BeamAlgorithmParams parameters;
	private int groupChannelMap; // 通道映射，表示通道分组的位图
	private int currentChanMap; // 当前通道映射，用于跟踪当前已处理的通道
	private int nextChanIndex = 0;
	private int firstSeqNum = 0;
	private int numChannels; // 通道数量，表示通道分组中的通道数
	private FFTDataUnit[] channelFFTUnits; // FFTDataUnit数组，用于存储通道分组的FFT数据单元
	private int arrayShape;
	private PamVector[] arrayMainAxes;
	
	public BeamGroupProcess(BeamFormerBaseProcess beamFormerBaseProcess, BeamAlgorithmProvider provider, 
			BeamAlgorithmParams parameters, int groupChannels, int firstSeqNum, int beamogramNum) {
		this.beamFormerBaseProcess = beamFormerBaseProcess;
		this.provider = provider;
		this.parameters = parameters;
		this.groupChannelMap = groupChannels;
		this.firstSeqNum = firstSeqNum;
		this.numChannels = PamUtils.getNumChannels(groupChannels);
		channelFFTUnits = new FFTDataUnit[numChannels];
		beamFormerAlgorithm = provider.makeAlgorithm(beamFormerBaseProcess, parameters, firstSeqNum, beamogramNum);
		beamFormerAlgorithm.prepare();
		

		// work out the subArray Shape. 
		ArrayManager arrayManager = ArrayManager.getArrayManager();
		PamArray currentArray = arrayManager.getCurrentArray();
		int phones = groupChannelMap;
		phones = beamFormerBaseProcess.getFftDataSource().getChannelListManager().channelIndexesToPhones(phones);
		arrayShape = arrayManager.getArrayShape(currentArray, phones);
		arrayMainAxes = arrayManager.getArrayDirections(currentArray, phones);
	}
	
	/**
	 * 
	 * @return The parameter set currently used for this set of channels. 
	 */
	public BeamAlgorithmParams getAlgorithmParams() {
		return parameters;
	}
	
	/**
	 * Process FFT data - will group up into correct channel groups
	 * @param fftDataUnit
	 */
	public void process(FFTDataUnit fftDataUnit) { // 读取数据存入到该分组的channel中
		int chMap = fftDataUnit.getChannelBitmap();
		channelFFTUnits[nextChanIndex++] = fftDataUnit;
		currentChanMap |= chMap;
		if (nextChanIndex == numChannels) { // 当索引等于该分组预设的通道数，算法才会开始处理
			beamFormerAlgorithm.process(channelFFTUnits);
			nextChanIndex = 0;
		}
	}

	/**
	 * @return the beamFormerAlgorithm
	 */
	public BeamFormerAlgorithm getBeamFormerAlgorithm() {
		return beamFormerAlgorithm;
	}

	/**
	 * @return the arrayShape
	 */
	public int getArrayShape() {
		return arrayShape;
	}

	/**
	 * @return the arrayMainAxes
	 */
	public PamVector[] getArrayMainAxes() {
		return arrayMainAxes;
	}

	/**
	 * Called from the localiser to reset the FFT store so 
	 * that it always starts from scratch as new data are sent
	 */
	public void resetFFTStore() {
		nextChanIndex = 0;
		currentChanMap = 0;
	}

	/**
	 * @return the groupChannelMap
	 */
	public int getGroupChannelMap() {
		return groupChannelMap;
	}

	/**
	 * @return the nextChanIndex
	 */
	public int getNextChanIndex() {
		return nextChanIndex;
	}

	/**
	 * @return the numChannels
	 */
	public int getNumChannels() {
		return numChannels;
	}

	/**
	 * @return the beamFormerBaseProcess
	 */
	public BeamFormerBaseProcess getBeamFormerBaseProcess() {
		return beamFormerBaseProcess;
	}
	
}

