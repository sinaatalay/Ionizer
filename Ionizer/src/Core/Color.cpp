#include "Ionpch.h"

#include "Color.h"
#include "Log.h"

namespace Ionizer {

	void Color::SetHSV(float fH, float fS, float fV) {
		if (fH > 360 || fH < 0 || fS>1 || fS < 0 || fV>1 || fV < 0) {
			LOG_CRITICAL("The given HSV values are not in the valid range!");
			return;
		}

		float fR, fG, fB;

		float fC = fV * fS; // Chroma
		float fHPrime = fmod(fH / 60.0, 6);
		float fX = fC * (1 - fabs(fmod(fHPrime, 2) - 1));
		float fM = fV - fC;

		if (0 <= fHPrime && fHPrime < 1) {
			fR = fC;
			fG = fX;
			fB = 0;
		}
		else if (1 <= fHPrime && fHPrime < 2) {
			fR = fX;
			fG = fC;
			fB = 0;
		}
		else if (2 <= fHPrime && fHPrime < 3) {
			fR = 0;
			fG = fC;
			fB = fX;
		}
		else if (3 <= fHPrime && fHPrime < 4) {
			fR = 0;
			fG = fX;
			fB = fC;
		}
		else if (4 <= fHPrime && fHPrime < 5) {
			fR = fX;
			fG = 0;
			fB = fC;
		}
		else if (5 <= fHPrime && fHPrime < 6) {
			fR = fC;
			fG = 0;
			fB = fX;
		}
		else {
			fR = 0;
			fG = 0;
			fB = 0;
		}

		fR += fM;
		fG += fM;
		fB += fM;

		r = (unsigned char)(fR * 255);
		g = (unsigned char)(fG * 255);
		b = (unsigned char)(fB * 255);
	}
	void Color::SetRGB(int R, int G, int B) {
		r = R;
		g = G;
		b = B;
	}

}

