package com.sanger.validateStructuralVariants;

import com.beust.jcommander.*;
import java.io.File;

public class FileValidator implements IParameterValidator {
	public void validate(String name, String value) throws ParameterException {
		File f = new File(value);
		if(!f.exists() || f.isDirectory())
			throw new ParameterException("Parameter " + name + " does not exist or is a directory. Found " + value);
	}
}