FluidGraphPlay : FluidRealTimeModel {
	var <>source, <>numBands, <>threshold, <>minDur, <>minDist,
    <>forget, <>start, <>output, <>windowSize, <>hopSize, <>fftSize, <>maxFFTSize;

		*new {|server, source = -1, numBands = 64, threshold = 0.3,
  minDur = 10, minDist = 10, forget = 1, start = 0, output,
		windowSize = 1024, hopSize = -1,
  fftSize = -1, maxFFTSize = 16384|
		^super.new(server,[source, numBands, threshold, minDur, minDist,
    forget, start, output, windowSize, hopSize, fftSize, maxFFTSize])
		.source_(source)
		.numBands_(numBands)
		.threshold_(threshold)
		.minDur_(minDur)
		.minDist_(minDist)
		.forget_(forget)
		.start_(start)
		.output_(output)
		.windowSize_(windowSize)
		.hopSize_(hopSize)
		.fftSize_(fftSize)
		.maxFFTSize_(maxFFTSize);
	}

	prGetParams{^[
		this.source, this.numBands,this.threshold, this.minDur, this.minDist,
		this.forget, this.start, this.output, this.windowSize,
		this.hopSize,this.fftSize, this.maxFFTSize,-1,-1];}

	analyze{|action|
		actions[\analyze] = [nil,action];
		this.prSendMsg(this.prMakeMsg(\analyze, id));
	}

	ar { arg start = 0, threshold = 0.1, minDur = 10, minDist = 10, forget = 100;
		source = source ?? {-1};
		output = output ?? {-1};
		^FluidGraphPlayQuery.ar(this, source, numBands, threshold, minDur, minDist,
    forget, start, output, windowSize, hopSize, fftSize, maxFFTSize);
	}

}



FluidGraphPlayQuery : UGen
{
	var <>pluginname;

	*ar { |...args|
        args = [1] ++ args.collect{|x| x.asUGenInput};
		^this.new1('audio',  "FluidGraphPlayQuery", *args)
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
