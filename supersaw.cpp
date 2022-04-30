#include <Bela.h>
#include <cmath>
#include <libraries/Trill/Trill.h>
#include <libraries/OnePole/OnePole.h>
#include <libraries/Oscillator/Oscillator.h>
#include <libraries/Biquad/Biquad.h>
#include <libraries/ADSR/ADSR.h>


#define NUM_TOUCH 5 // Number of touches on Trill sensor

// filter stuff
Biquad lpFilter[NUM_TOUCH];	// Biquad low-pass frequency;
Biquad hpFilter;	// Biquad high-pass frequency;

float gHPfreq = 80.0;	// Cut-off frequency for high-pass filter (Hz)
float gLPfreq = 200.0;	// Cut-off frequency for low-pass filter (Hz)

float gFilterQ = 0.707; // Quality factor for the biquad filters to provide a Butterworth response

// setup envelopes
ADSR envelope[NUM_TOUCH]; // ADSR envelope

float gAttack = 1.1; // Envelope attack (seconds)
float gDecay = 0.25; // Envelope decay (seconds)
float gRelease = 5.0; // Envelope release (seconds)
float gSustain = 1.0; // Envelope sustain level


// Trill object declaration
Trill touchSensor;

// Location of touches on Trill Bar
float gTouchLocation[NUM_TOUCH] = { 0.0, 0.0, 0.0, 0.0, 0.0 };
// Size of touches on Trill Bar
float gTouchSize[NUM_TOUCH] = { 0.0, 0.0, 0.0, 0.0, 0.0 };
// Number of active touches
int gNumActiveTouches = 0;

// Sleep time for auxiliary task in microseconds
unsigned int gTaskSleepTime = 12000; // microseconds

// One Pole filters objects declaration
OnePole freqFilt[NUM_TOUCH], ampFilt[NUM_TOUCH];
// Frequency of one pole filters
float gCutOffFreq = 5, gCutOffAmp = 15;
// Oscillators objects declaration
Oscillator osc[NUM_TOUCH];
// Range for oscillator frequency mapping
float gFreqRange[2] = { 88.0, 92.0 };
// Range for oscillator amplitude mapping
float gAmplitudeRange[2] = { 0.0, 1.0 } ;

/*
* Function to be run on an auxiliary task that reads data from the Trill sensor.
* Here, a loop is defined so that the task runs recurrently for as long as the
* audio thread is running.
*/
void loop(void*)
{
	while(!Bela_stopRequested())
	{
		// Read locations from Trill sensor
		touchSensor.readI2C();
		gNumActiveTouches = touchSensor.getNumTouches();
		for(unsigned int i = 0; i <  gNumActiveTouches; i++) {
			gTouchLocation[i] = touchSensor.touchLocation(i);
			gTouchSize[i] = touchSensor.touchSize(i);
		}
		// For all inactive touches, set location and size to 0
		for(unsigned int i = gNumActiveTouches; i < NUM_TOUCH; i++) {
			gTouchLocation[i] = 0.0;
			gTouchSize[i] = 0.0;
		}
		usleep(gTaskSleepTime);
	}
}

bool setup(BelaContext *context, void *userData)
{
	// Setup a Trill Bar sensor on i2c bus 1, using the default mode and address
	if(touchSensor.setup(1, Trill::BAR) != 0) {
		fprintf(stderr, "Unable to initialise Trill Bar\n");
		return false;
	}
	touchSensor.printDetails();

	// Set and schedule auxiliary task for reading sensor data from the I2C bus
	Bela_runAuxiliaryTask(loop);

	// For each possible touch...
	for(unsigned int i = 0; i < NUM_TOUCH; i++) {
		// Setup corresponding oscillator
		osc[i].setup(context->audioSampleRate, Oscillator::sawtooth);
		// Setup low pass filters for smoothing frequency and amplitude
		freqFilt[i].setup(gCutOffFreq, context->audioSampleRate);
		ampFilt[i].setup(gCutOffAmp, context->audioSampleRate);
		// setup envelope
		envelope[i].setAttackRate(gAttack * context->audioSampleRate);
		envelope[i].setDecayRate(gDecay * context->audioSampleRate);
		envelope[i].setReleaseRate(gRelease * context->audioSampleRate);
		envelope[i].setSustainLevel(gSustain);
	}

	// setup filter
	Biquad::Settings settings{
			.fs = context->audioSampleRate,
			.cutoff = gLPfreq,
			.type = Biquad::lowpass,
			.q = gFilterQ,
			.peakGainDb = 0,
			};
	for (unsigned int i = 0; i < NUM_TOUCH; i++) {
		lpFilter[i].setup(settings);
	}
	settings.cutoff = gHPfreq;
	settings.type = Biquad::highpass;
	hpFilter.setup(settings);



	return true;
}

void render(BelaContext *context, void *userData)
{
	for(unsigned int i = 0; i < NUM_TOUCH; i++) {
		if (gTouchSize[i]>0) {
			envelope[i].gate(true);
		} else {
			envelope[i].gate(false);
		}
		lpFilter[i].setFc(freqFilt[i].process(2.0*pow(10,map(gTouchLocation[i], 0, 1, 1,4))));
	}
	for(unsigned int n = 0; n < context->audioFrames; n++) {
		float out = 0.0;
		/* For each touch:
		*
		* 	- Map touch location to frequency of the oscillator
		* 	and smooth value changes using a single pole LP filter
		* 	- Map touch size toa amplitude of the oscillator and
		* 	smooth value changes using a single pole LP filter
		* 	- Compute oscillator value and add to output.
		* 	- The overall output will be scaled by the number of touches.
		*/
		for(unsigned int i = 0; i < NUM_TOUCH; i++) {
			float frequency, amplitude;
			// if (gTouchLocation[i]<0.25) {
			// 	frequency = map(gTouchLocation[i], 0, 0.25, gFreqRange[0], gFreqRange[1])*2;
			// } else if (gTouchLocation[i]<0.5) {
			// 	frequency = map(gTouchLocation[i], 0.25, 0.5, gFreqRange[0], gFreqRange[1])/2;
			// } else if (gTouchLocation[i]<0.75) {
			// 	frequency = map(gTouchLocation[i], 0.5, 0.75, gFreqRange[0], gFreqRange[1]);
			// } else {
			// 	frequency = map(gTouchLocation[i], 0.75, 1, gFreqRange[0], gFreqRange[1])*4;
			// }
			// frequency = freqFilt[i].process(frequency);
			// amplitude = map(gTouchSize[i], 0, 1, gAmplitudeRange[0], gAmplitudeRange[1]);
			frequency = 90;
			if (i>0) {
				frequency=frequency+(gTouchLocation[i]-gTouchLocation[i-1])*10;
			}
			amplitude = envelope[i].process()/NUM_TOUCH;

			out += lpFilter[i].process((1.f/NUM_TOUCH) * amplitude * osc[i].process(frequency));
		}

		// Process input signal with high pass filter
		out = hpFilter.process(out);

		// Write computed output to audio channels
		for(unsigned int channel = 0; channel < context->audioOutChannels; channel++) {
			audioWrite(context, n, channel, out);
		}
	}
}

void cleanup(BelaContext *context, void *userData)
{}
