
typedef uint64_t truth6;

struct NPN {
	bool oc = false;
	bool ic[6] = {};
	int p[6] = {-1, -1, -1, -1, -1, -1};

	static NPN identity(int ninputs) {
		NPN ret;
		for (int i = 0; i < ninputs; i++)
			ret.p[i] = i;
		return ret;
	}

	bool is_identity() const
	{
		if (oc)
			return false;
		for (int i = 0; i < 6 && p[i] != -1; i++)
		if (p[i] != i || ic[i])
			return false;
		return true;
	}

	int ninputs() const
	{
		int j;
		for (j = 0; j < 6 && p[j] != -1; j++);
		return j;
	}

	truth6 operator()(truth6 m) const {
		truth6 ret = 0;
		for (int idx1 = 0; idx1 < (1 << ninputs()); idx1++) {
			if (!(m & (truth6) 1 << idx1))
				continue;
			int idx2 = 0;
			for (int j = 0; j < ninputs(); j++) {
				if ((!!(idx1 & 1 << j)) ^ ic[j])
					idx2 |= 1 << p[j];
			}
			ret |= (truth6) 1 << idx2;
		}
		return oc ? ~ret : ret;
	}

	NPN operator*(const NPN &other)
	{
		NPN ret;
		ret.oc = oc ^ other.oc;
		for (int i = 0; i < 6 && p[i] != -1; i++) {
			ret.p[i] = p[other.p[i]];
			ret.ic[i] = ic[other.p[i]] ^ other.ic[i];
		}
		return ret;
	}

	void dump()
	{
		printf("%s(", oc ? "#" : "");
		for (int i = 0; i < 6 && p[i] != -1; i++)
			printf("%s%d->%d%s", i != 0 ? " " : "", i, p[i], ic[i] ? "#" : "");
		printf(")\n");
	}

	NPN inv()
	{
		NPN ret;
		ret.oc = oc;
		for (int i = 0; i < 6 && p[i] != -1; i++) {
			ret.p[p[i]] = i;
			ret.ic[p[i]] = ic[i];
		}
		return ret;
	}

	int c_fingerprint()
	{
		return oc | ic[0] << 1 | ic[1] << 2 | ic[2] << 3 | ic[3] << 4 | ic[4] << 5 | ic[5] << 6;
	}
};

truth6 npn_semiclass(truth6 m, int ninputs, NPN &npn);
void npn_semiclass_allrepr(truth6 m, int ninputs,
						   std::function<void(truth6, NPN&)> cb);
