import deepspeech
import numpy as np
from scipy.io.wavfile import read

# Load the DeepSpeech model
model = deepspeech.Model("deepspeech-0.9.3-models.pbmm")

def transcribe_audio(file_path):
    fs, audio = read(file_path)
    audio = np.frombuffer(audio, dtype=np.int16)
    transcript = model.stt(audio)
    print("Transcript:", transcript)

transcribe_audio("path_to_audio_file.wav")

