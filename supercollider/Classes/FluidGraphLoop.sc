FluidGraphLoop : FluidRealTimeModel {
	var <>source, <>numBands, <>threshold,
	<>quantize, <>start, <>end,
	<>output, <>windowSize, <>hopSize, <>fftSize, <>maxFFTSize;

		*new {|server, source = -1, numBands = 64, threshold = 0.3,
  quantize = 0, start = 0, end = 1, output, windowSize = 1024, hopSize = -1,
  fftSize = -1, maxFFTSize = 16384|
		^super.new(server,[source, numBands, threshold, quantize, start,
    end, output,windowSize, hopSize, fftSize, maxFFTSize])
		.source_(source)
		.numBands_(numBands)
		.threshold_(threshold)
		.quantize_(quantize)
		.start_(start)
		.end_(end)
		.output_(output)
		.windowSize_(windowSize)
		.hopSize_(hopSize)
		.fftSize_(fftSize)
		.maxFFTSize_(maxFFTSize);
	}

	prGetParams{^[
		this.source, this.numBands,this.threshold, this.quantize, this.start,
		this.end, this.output, this.windowSize,
		this.hopSize,this.fftSize, this.maxFFTSize,-1,-1];}

	analyze{|action|
		actions[\analyze] = [nil,action];
		this.prSendMsg(this.prMakeMsg(\analyze, id));
	}

	ar { arg start = 0, end = 1;
		source = source ?? {-1};
		output = output ?? {-1};
		^FluidGraphLoopQuery.ar(this, source, numBands, threshold, quantize, start,
    end, output,windowSize, hopSize, fftSize, maxFFTSize);
	}

}



FluidGraphLoopQuery : UGen
{
	var <>pluginname;

	*ar { |...args|
        args = [1] ++ args.collect{|x| x.asUGenInput};
		^this.new1('audio',  "FluidGraphLoopQuery", *args)
	}

	init { |pluginname...args|
		this.pluginname = pluginname;
		inputs = args;
		rate = 'audio';
		specialIndex = 0;
	}


	name{
		^pluginname.asString;
	}
}
