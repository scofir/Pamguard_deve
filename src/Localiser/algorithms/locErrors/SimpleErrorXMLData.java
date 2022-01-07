package Localiser.algorithms.locErrors;

import java.text.DecimalFormat;

public class SimpleErrorXMLData extends ErrorXMLData {

	private String errorType;
	private double[] errorData;
	
	
	public SimpleErrorXMLData(String errorType, double[] errorData) {
		this.errorType = errorType;
		this.errorData = errorData;
	}
	
	/* (non-Javadoc)
	 * @see Localiser.algorithms.locErrors.ErrorXMLData#getXMLString()
	 */
	@Override
	public String getXMLString() {
		String x = String.format("{\"Name\":\"%s\",\"Data\":[", errorType);
		for (int i = 0; i < errorData.length-1; i++) {
			x += errorDigitFormat.format(errorData[i]) + ",";
		}
		x += errorDigitFormat.format(errorData[errorData.length-1]) + "]}";
		return x;
	}

}
