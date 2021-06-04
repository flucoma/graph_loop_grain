FluidGraphGrain : FluidRealTimeModel {
	var <>source, <>numBands, <>threshold,
	<>numClusters, <>forgetfulness, <>randomness, <>phase, <>start,
	<>output, <>windowSize, <>hopSize, <>fftSize, <>maxFFTSize;

		*new {|server, source = -1, numBands = 64, threshold = 0.3,
  numClusters = 10, forgetfulness = 100, randomness = 0.1,
  phase = 1, start = 0, output, windowSize = 1024, hopSize = -1,
  fftSize = -1, maxFFTSize = 16384|
		^super.new(server,[source, numBands, threshold, numClusters, forgetfulness,
    randomness, phase, start, output,windowSize, hopSize, fftSize, maxFFTSize])
		.source_(source)
		.numBands_(numBands)
		.threshold_(threshold)
		.numClusters_(numClusters)
		.forgetfulness_(forgetfulness)
		.randomness_(randomness)
		.phase_(phase)
		.start_(start)
		.output_(output)
		.windowSize_(windowSize)
		.hopSize_(hopSize)
		.fftSize_(fftSize)
		.maxFFTSize_(maxFFTSize);
	}

	prGetParams{^[
		this.source, this.numBands,this.threshold, this.numClusters, this.forgetfulness,
		this.randomness, this.phase, this.start, this.output, this.windowSize,
		this.hopSize,this.fftSize, this.maxFFTSize,-1,-1];}

	analyze{|action|
		actions[\analyze] = [nil,action];
		this.prSendMsg(this.prMakeMsg(\analyze, id));
	}

	ar { arg start = 0, threshold = 0.1, forgetfulness = 100, randomness = 0.1, phase = 1;
		source = source ?? {-1};
		output = output ?? {-1};
		^FluidGraphGrainQuery.ar(this, source, numBands, threshold, numClusters, forgetfulness,
			randomness, phase, start, output, windowSize, hopSize, fftSize, maxFFTSize);
	}

}



FluidGraphGrainQuery : UGen
{
	var <>pluginname;

	*ar { |...args|
        args = [1] ++ args.collect{|x| x.asUGenInput};
		^this.new1('audio',  "FluidGraphGrainQuery", *args)
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
