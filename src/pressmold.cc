//
// pressmold -- an OpenSTA-based standard cell mapper
//

#define CUT_MAXIMUM		6

//#define SIBLING_RECORDING

#include <sta/Sta.hh>
#include <sta/StaMain.hh>
#include <sta/Network.hh>
#include <sta/Liberty.hh>
#include <sta/FuncExpr.hh>
#include <sta/PortDirection.hh>
#include <sta/ConcreteNetwork.hh>
#include <sta/VerilogNamespace.hh>
#include <sta/TimingArc.hh>
#include <sta/TimingRole.hh>
#include <sta/DcalcAnalysisPt.hh>
#include <sta/ArcDelayCalc.hh>
#include <sta/Corner.hh>

#include <tcl.h>
#include <tclreadline.h>

#include <algorithm>
#include <iostream>
#include <fstream>
#include <random>
#include <vector>
#include <map>
#include <cstdlib>
#include <cstdint>
#include <cassert>

#include "npn.h"

struct TargetIndex {
	struct {
		sta::LibertyCell *cell;
		sta::LibertyPort *hi;
		sta::LibertyPort *lo;
	} tie;

	sta::LibertyCell *inv_cell;

	struct Target {
		sta::LibertyCell *cell;
		NPN map;
	};

	std::map<std::pair<truth6, int>, std::vector<Target>> classes;
} target_index;

using Target = TargetIndex::Target;

struct AndNode;
struct CutList {
	AndNode **array;
	int size;

	CutList(AndNode **array)
		: array(array)
	{
		AndNode **p;
		for (p = array; (p < array + CUT_MAXIMUM) && (*p != NULL); p++);
		size = p - array;
	}

	CutList(std::vector<AndNode *> &vector)
		: array(&vector[0]), size(vector.size())
	{
	}

	class iterator {
		AndNode **p;
	public:
		typedef std::input_iterator_tag iterator_category;
		typedef AndNode* value_type;
		typedef ptrdiff_t difference_type;
		typedef AndNode** pointer;
		typedef AndNode*& reference;

		iterator(AndNode **p) : p(p) {}
		iterator& operator++() { p++; return *this; }
		bool operator==(const iterator &other) const
			{ return p == other.p; }
		bool operator==(AndNode **p) const
			{ return this->p == p; }
		int operator-(const iterator &other) const
			{ return p - other.p; }
		bool operator!=(const iterator &other) const { return !(*this == other); }
		AndNode *operator*() const { return *p; }
	};
	iterator begin() const { return iterator(array); };
	iterator end()   const { return iterator(array + size); };
};

struct NodeInput {
	NodeInput() {}
	NodeInput(AndNode *node, bool negated)
	{
		set_node(node);
		this->negated = negated;
	}

	AndNode *node = NULL;
	bool negated = false;

	void negate()					{ negated ^= true; }
	void set_node(AndNode *source)	{ negated = false; node = source; }

	void set_const(int state)
	{
		set_node(NULL);
		assert(state == 1 || state == 0);
		if (state)
			negate();
	}

	bool is_const() { return !node; }
	bool eval()		{ assert(is_const()); return negated; }
	bool polarity();

	bool operator<(const NodeInput &other) const
	{
		if (!node != !other.node)
			return !node < !other.node;
		return std::tie(negated, node) < std::tie(other.negated, other.node);
	}
	bool operator==(const NodeInput &other) const
			{ return std::tie(negated, node) == std::tie(other.negated, other.node); }

	unsigned int hash() const
	{
		return (uintptr_t) node + negated;
	}

	std::vector<bool> truth_table(std::vector<AndNode *> &cut);
};

struct AndNode {
	bool pi = false;
	bool po = false;

	NodeInput ins[2];
	AndNode *sibling = NULL;
	bool polarity = false;
	bool in_repr = false;

	std::string label;
	AndNode() {};

	struct Match {
		truth6 semiclass;
		NPN npn;
		AndNode *cut[CUT_MAXIMUM];
#ifdef SIBLING_RECORDING
		std::vector<AndNode *> used_siblings;
#endif
	};

	// Scratch area for algorithms
	bool visited;
	union {
		int idx;
		int refs;
		struct {
			Match *matches;
			struct {
				union {
					float farea;
					int depth;
				};
				float flow_fouts;
				float area;
				float fuzzy_fouts;

				// an invariant: if matches_valid on Network is true
				// and this node's map_fouts is non-zero, then sel/sel_target
				// must be valid
				int sel;
				Target *sel_target;

				int map_fouts;
				int save_sel;
				Target *save_sel_target;
				float save_area;
				sta::Net *net;
			} pol[2];
		};
		AndNode *replacement;
	};

	sta::Net *net;
	int fanouts;
	int fid; // frontier index

	void apply_replacements()
	{
		for (int i = 0; i < 2; i++)
		if (ins[i].node && ins[i].node->replacement) {
			assert(ins[i].node != ins[i].node->replacement);
			ins[i].node = ins[i].node->replacement;
		}

		if (sibling && sibling->replacement)
			sibling = sibling->replacement;
	}

	bool expand();

	struct FaninList {
		AndNode *node;
		bool include_sibling;
		class iterator {
			AndNode *node; int pos;
		public:
			typedef std::input_iterator_tag iterator_category;
			typedef AndNode* value_type;
			typedef ptrdiff_t difference_type;
			typedef AndNode** pointer;
			typedef AndNode*& reference;

			iterator(AndNode *node, int pos)
				: node(node), pos(pos) {}
			iterator& operator++()
				{ pos++; return normalize(); }
			iterator& normalize()
				{ while (pos < 2 && !*(*this)) pos++; return *this; }
			bool operator==(const iterator &other) const
				{ return std::tie(node, pos) == std::tie(other.node, other.pos); }
			bool operator!=(const iterator &other) const
				{ return !(*this == other); }
			AndNode *operator*() const
				{ return pos == -1 ? node->sibling : node->ins[pos].node; }
		};
		iterator begin() const
			{ return iterator(node, include_sibling ? -1 : 0).normalize(); };
		iterator end()   const
			{ return iterator(node, 2); };
	};

	FaninList fanins()		{ return FaninList{ this, false }; }
	FaninList pointees()	{ return FaninList{ this, true }; }

	bool detect_mux(AndNode* &s, AndNode* &a, AndNode* &b);

	std::vector<bool> truth_table(std::vector<AndNode *> &cut, bool negate=false)
	{
		int index = 0;
		for (auto cut_node : cut) {
			if (cut_node == this) {
				std::vector<bool> ret;
				for (int i = 0; i < (1 << cut.size()); i++)
					ret.push_back((((i >> index) & 1) ^ negate) != 0);
				return ret;
			}
			index++;
		}

		if (pi)
			return {};

		std::vector<bool> a = ins[0].truth_table(cut);
		std::vector<bool> b = ins[1].truth_table(cut);

		if (a.empty() || b.empty())
			return {};

		assert(a.size() == 1 << cut.size() && b.size() == 1 << cut.size());
		std::vector<bool> ret;
		for (int i = 0; i < (1 << cut.size()); i++)
			ret.push_back((a[i] && b[i]) ^ negate);
		return ret;
	}
};

std::vector<bool> NodeInput::truth_table(std::vector<AndNode *> &cut)
{
	if (!node) {
		return std::vector<bool>(1 << cut.size(), negated);
	} else {
		return node->truth_table(cut, negated);
	}
}

bool NodeInput::polarity()
{
	if (is_const())
		return eval();
	else
		return node->polarity ^ negated;
}

bool AndNode::detect_mux(AndNode* &s, AndNode* &a, AndNode* &b)
{
	AndNode *n1 = ins[0].node, *n2 = ins[1].node;
	if (!n1 || !n2)
		return false;
	if (!n1->ins[0].node || !n1->ins[1].node)
		return false;
	if (!n2->ins[0].node || !n2->ins[1].node)
		return false;

	// TODO: has some rare false positives
	if (n1->ins[0].node == n2->ins[0].node) {
		s = n1->ins[0].node;
		a = n1->ins[1].node;
		b = n2->ins[1].node;
	} else if (n1->ins[0].node == n2->ins[1].node) {
		s = n1->ins[0].node;
		a = n1->ins[1].node;
		b = n2->ins[0].node;
	} else if (n1->ins[1].node == n2->ins[0].node) {
		s = n1->ins[1].node;
		a = n1->ins[0].node;
		b = n2->ins[1].node;
	} else if (n1->ins[1].node == n2->ins[1].node) {
		s = n1->ins[1].node;
		a = n1->ins[0].node;
		b = n2->ins[0].node;
	} else {
		return false;
	}
	return true;
}

truth6 mask6(int n) {
	if (n == 6)
		return ~(truth6) 0;
	return (((truth6) 1) << (1 << n)) - 1;
}

truth6 recode6(truth6 t1, CutList vars1, CutList vars2) {
	truth6 t2 = 0;
	for (int i2 = 0; i2 < 1 << vars2.size; i2++) {
		int i1 = 0;
		auto var2it = vars2.begin();
		int n1 = 0;
		for (auto var1 : vars1) {
			while (*var2it < var1)
				++var2it;
			assert(var2it != vars2.end() && *var2it == var1);
			if (i2 & 1 << (var2it - vars2.begin()))
				i1 |= 1 << n1;
			n1++;
		}
		if (t1 & (truth6) 1 << i1)
			t2 |= (truth6) 1 << i2;
	}
	return t2;
}

truth6 reduce6(truth6 t1, int nvars, uint32_t &removed)
{
	removed = (1 << nvars) - 1;

	for (int i = 0; i < 1 << nvars; i++) {
		for (int j = 0; j < nvars; j++) {
			if (i & (1 << j))
				continue;
			if ((!!(t1 & (truth6) 1 << i)) != !!(t1 & (truth6) 1 << (i | 1 << j)))
				removed &= ~(1 << j);
		}
	}

	if (!removed)
		return t1;

	truth6 t2 = 0;
	for (unsigned int i1 = 0, i2 = 0; i1 < (1 << nvars);
			i1 = ((i1 | removed) + 1) & ~removed, i2++) {
		if (t1 & (truth6) 1 << i1)
			t2 |= (truth6) 1 << i2;
	}

	return t2;
}

uint32_t read_be32(std::istream &f) {
	return ((uint32_t) f.get() << 24) |
		((uint32_t) f.get() << 16) | 
		((uint32_t) f.get() << 8) | (uint32_t) f.get();
}

struct Network {
	std::string name;
	std::vector<AndNode> node_storage;
	std::vector<std::unique_ptr<AndNode::Match[]>> match_storage;
	bool matches_valid = false;

	// We want to be code-compatible with toymap in how we iterate
	// over nodes. For that reason we have the NodeList indirection
	// which lets us type
	//
	//   for (auto node : nodes)
	//     node->...
	//
	struct NodeList {
		std::vector<AndNode> &storage;
		bool indication = false;

		class iterator {
			const NodeList &ms;
			AndNode *node;
		public:
			typedef std::input_iterator_tag iterator_category;
			typedef AndNode* value_type;
			typedef ptrdiff_t difference_type;
			typedef AndNode** pointer;
			typedef AndNode*& reference;

			iterator(const NodeList &ms, AndNode *node)
				: ms(ms), node(node) {}
			iterator& operator++() { node++; update(); return *this; }
			bool operator==(const iterator &other) const
				{ return node == other.node; }
			bool operator!=(const iterator &other) const
				{ return !(*this == other); }
			AndNode* operator*() const
				{ return node; }
			void update() {
				if (ms.indication) {
					int no = node - &ms.storage.front();
					if (no == ms.storage.size()) {
						printf("                  \r");
						fflush(stdout);
					} else if (!(no % 100)) {
						printf(" %6d/%6d\r", no, (int) ms.storage.size());
						fflush(stdout);
					}
				}
			}
		};
		iterator begin() const
			{ return iterator(*this, &storage.front()); };
		iterator end()   const
			{ return iterator(*this, &storage.front() + storage.size()); };
		class riterator {
			AndNode *node;
		public:
			typedef std::input_iterator_tag iterator_category;
			typedef AndNode* value_type;
			typedef ptrdiff_t difference_type;
			typedef AndNode** pointer;
			typedef AndNode*& reference;

			riterator(AndNode *node)
				: node(node) {}
			riterator& operator++() { node--; return *this; }
			bool operator==(const riterator &other) const
				{ return node == other.node; }
			bool operator!=(const riterator &other) const
				{ return !(*this == other); }
			AndNode* operator*() const
				{ return node; }
		};
		riterator rbegin() const
			{ return riterator(&storage.back()); };
		riterator rend()   const
			{ return riterator(&storage.front() - 1); };
		AndNode* operator[](int index)
			{ return &storage[index]; }
		int size() const
			{ return storage.size(); }
		NodeList w_indication()
			{ return NodeList{storage, true}; }
	} nodes;

	Network() : nodes(node_storage) {}
	~Network() {}

	Network (const Network&) = delete;
	Network& operator= (const Network&) = delete;
	Network(Network&& other)
			: nodes(node_storage) {
		name = other.name;
		std::swap(node_storage, other.node_storage);
		std::swap(match_storage, other.match_storage);
		matches_valid = other.matches_valid;
	}
	Network& operator=(Network&& other) {
		name = other.name;
		std::swap(node_storage, other.node_storage);
		std::swap(match_storage, other.match_storage);
		matches_valid = other.matches_valid;
		return *this;
	}

	static Network read_aiger(std::ifstream &f, sta::ConcreteNetwork *stan=NULL,
							  const char *name="top", const char *filename="")
	{
		assert(f.get() == 'a' && f.get() == 'i'
			   && f.get() == 'g' && f.get() == ' ');

		int M, I, L, O, A;
		f >> M >> I >> L >> O >> A;
		assert(!f.fail() && !f.eof());
		assert(f.get() == '\n');
		assert(L == 0); // no latches

		Network ret;
		ret.name = name;
		ret.node_storage.resize(I + A + O);
		auto &nodes = ret.nodes;

		for (int j = 0; j < I; j++)
			nodes[j]->pi = true;

		for (int j = I + A; j < I + A + O; j++) {
			AndNode *node = nodes[j];
			node->po = true;
			node->ins[1].set_const(1);

			int pivot;
			f >> pivot;
			assert(f.get() == '\n');
			assert(pivot >= 0 && pivot <= (int) nodes.size() * 2);
			NodeInput &in = node->ins[0];
			if (pivot < 2) {
				in.set_const(pivot);
			} else {
				in.set_node(nodes[pivot / 2 - 1]);
				if (pivot & 1)
					in.negate();
			}
		}

		for (int j = 0; j < A; j++) {
			int pivot = 2 * (I + j) + 2;
			AndNode *node = nodes[I + j];

			for (int p = 0; p < 2; p++) {
				int c, delta = 0, shift = 0;
				while ((c = f.get()) != EOF) {
					delta |= (c & 0x7f) << shift;
					shift += 7;
					if (!(c & 0x80))
						break;
				}
				pivot = pivot - delta;
				assert(pivot >= 0 && pivot <= (int) nodes.size() * 2);
				NodeInput &in = node->ins[p];
				if (pivot < 2) {
					in.set_const(pivot);
				} else {
					in.set_node(nodes[pivot / 2 - 1]);
					if (pivot & 1)
						in.negate();
				}
			}
			node->polarity = node->ins[0].polarity() &&
								node->ins[1].polarity();
		}

		// update POs polarities once the guts were filled in
		for (int j = I + A; j < I + A + O; j++) {
			AndNode *node = nodes[j];
			node->polarity = node->ins[0].polarity() &&
								node->ins[1].polarity();
		}

		int c;
		while ((c = f.get()) != EOF) {
			if (c == 'i' || c == 'o') {
				int i; std::string s;
				f >> i >> s;
				assert(!f.eof() && !f.fail());
				if (c == 'o')
					i += I + A;
				assert(i >= 0 && i < (int) nodes.size());
				nodes[i]->label = "\\" + s;
			} else if (c == 'c') {
				break;
			} else if (c == '\n') {
				continue;
			} else {
				printf("unhandled '%c'\n", c);
			}
		}

		int ni = 0, no = 0;
		for (auto node : nodes) {
			if (node->pi && node->label.empty()) {
				char name[30];
				snprintf(name, sizeof(name), "\\i%04d", ni++);
				node->label = name;
			}
			if (node->po && node->label.empty()) {
				char name[30];
				snprintf(name, sizeof(name), "\\o%04d", no++);
				node->label = name;
			}

			node->sibling = NULL;
		}

		if (ni || no)
			printf("Made up %d input and %d output names\n", ni, no);

		while ((c = f.get()) != EOF) {
			switch (c) {
			case '\n':
				goto done;
			case 'q':
				{
					read_be32(f);
					int pairnum = read_be32(f);
					for (int i = 0; i < pairnum; i++) {
						int repr = read_be32(f);
						int sibling = read_be32(f);
						assert(repr >= 0 && sibling >= 0);
						assert(repr - 1 < nodes.size() && sibling < repr);
						assert(!f.eof() && !f.fail());
						if (sibling == 0) {
							printf("Warning: constant choice! Ignoring.\n");
							continue;
						}
						nodes[repr - 1]->sibling = nodes[sibling - 1];
					}
				}
				break;
			case 'n':
				{
					uint32_t len = read_be32(f);
					ret.name = std::string(len, '\0');
					f.read(&ret.name[0], len);

					// TODO: we are resetting the name because it usually
					// contains slashes
					ret.name = name;
				}
			default:
				{
					uint32_t len = read_be32(f);
					f.ignore(len);
					printf("section '%c' (%d): ignoring\n", c, c); // %d bytes\n", c, c, len);
				}
				break;
			}
		}
		done:

		while (true) {
			std::string scratch;
			std::getline(f, scratch);
			if (f.eof())
				break;
			assert(!f.fail());
			printf("input file: %s\n", scratch.c_str());
		}

		ret.verify();
		printf("Read network '%s' with %d nodes\n", ret.name.c_str(), A);

		if (stan) {
			sta::Library *lib = stan->findLibrary("mapping");
			if (!lib)
				lib = stan->makeLibrary("mapping", NULL);
			sta::Cell *network_cell = stan->findCell(lib, ret.name.c_str());
			if (network_cell)
				stan->deleteCell(network_cell);
			network_cell = stan->makeCell(lib, ret.name.c_str(), false, filename);

			for (auto node : ret.nodes)
			if (node->pi || node->po) {
				std::string port_name = sta::portVerilogToSta(node->label.c_str());
				sta::Port *port = stan->makePort(network_cell, port_name.c_str());
				stan->setDirection(port, node->pi ? sta::PortDirection::input()
											: sta::PortDirection::output());
			}
		}

		return ret;
	}

	void fanouts(bool discount_mux=false, bool repr_usage=true)
	{
		for (auto node : nodes) {
			node->in_repr = false;
			node->fanouts = 0;
		}

		if (repr_usage)
		for (auto it = nodes.rbegin(); it != nodes.rend(); ++it) {
			if ((*it)->in_repr || (*it)->po)
			for (auto fanin : (*it)->fanins())
				fanin->in_repr = true;
		}

		for (auto node : nodes) {
			if (node->po)
				node->fanouts++;

			if (repr_usage && !node->in_repr)
				continue;

			for (auto fanin : node->fanins())
				fanin->fanouts++;

			AndNode *s, *a, *b;
			if (discount_mux && node->detect_mux(s, a, b)) {
				assert(s->fanouts-- >= 1);
				if (a == b)
					assert(a->fanouts-- >= 1);
			}
		}
	}

	void verify()
	{
		fanouts();

		for (auto node : nodes) {
			if (node->sibling)
				assert(!node->sibling->fanouts);

			if (node->pi)
				assert(!node->polarity);
			else
				assert(node->polarity == (node->ins[0].polarity() && node->ins[1].polarity()));
		}
	}

	void stats()
	{
		int nand_nodes = 0, nins = 0, nouts = 0, nmuxes = 0, nsiblings = 0;
		for (auto node : nodes) {
			if (node->pi) {
				nins++;
			} else if (node->po) {
				nouts++;
			} else {
				nand_nodes++;
				AndNode *a, *b, *s;
				if (node->detect_mux(s, a, b))
					nmuxes++;
			}

			if (node->sibling)
				nsiblings++;
		}
		printf("Mapping problem summary:\n");
		printf("  %d inputs %d outputs %d nodes %d choice pairs %d xor/mux detections\n",
				nins, nouts, nand_nodes, nsiblings, nmuxes);
	}

	void dump_verilog(const char *module_name, std::ostream &f)
	{
		f << "module \\" << module_name << " (";
		bool first = true;
		for (auto node : nodes)
		if (node->pi || node->po) {
			if (first)
				first = false;
			else 
				f << ", ";
			f << node->label.c_str() << " ";
		}
		f << ");\n";
		int idx = 0;
		char scratch[64];
		for (auto node : nodes) {
			node->idx = idx++;
			if (node->pi) {
				snprintf(scratch, sizeof(scratch), "  input wire %s ;\n",
						 node->label.c_str());
				f << scratch;
				snprintf(scratch, sizeof(scratch), "  wire $%08d = %s ;\n",
						 node->idx, node->label.c_str());
				f << scratch;
				continue;
			}
			if (node->ins[0].node && node->ins[1].node) {
				snprintf(scratch, sizeof(scratch), "  wire $%08d = %s$%08d && %s$%08d;\n",
						 node->idx, node->ins[0].negated ? "!" : "", node->ins[0].node->idx,
						 node->ins[1].negated ? "!" : "", node->ins[1].node->idx);
			} else if (node->ins[0].node && !node->ins[1].node) {
				assert(node->ins[1].eval());
				snprintf(scratch, sizeof(scratch), "  wire $%08d = %s$%08d;\n",
						 node->idx, node->ins[0].negated ? "!" : "", node->ins[0].node->idx);
			} else if (!node->ins[0].node && !node->ins[1].node) {
				snprintf(scratch, sizeof(scratch), "  wire $%08d = %d;\n",
						 node->idx, node->ins[0].eval() && node->ins[1].eval());
			} else {
				abort();
			}
			f << scratch;
			if (node->po) {
				snprintf(scratch, sizeof(scratch), "  output wire %s = $%08d;\n",
						 node->label.c_str(), node->idx);
				f << scratch;
			}
		}
		f << "endmodule\n";
	}

	int consolidate()
	{
		// pointers within cuts are invalidated by the move
		invalidate_matches();

		std::vector<AndNode*> used;

		for (auto node : nodes) {
			node->replacement = NULL;
			node->refs = 0;

			if (node->po) {
				assert(node->refs++ == 0);
				used.push_back(node);
			}
		}

		for (int i = 0; i < (int) used.size(); i++) {
			for (auto pointee : used[i]->pointees()) {
				if (!pointee->refs++)
					used.push_back(pointee);
				assert(pointee->refs > 0);
			}
		}

		used.clear();

		for (auto node : nodes) {
			if (node->po) {
				assert(--node->refs == 0);
				used.push_back(node);
			}
		}

		for (int i = 0; i < (int) used.size(); i++) {
			for (auto pointee : used[i]->pointees()) {
				if (!--pointee->refs)
					used.push_back(pointee);
				assert(pointee->refs >= 0);
			}
		}

		for (auto it = used.rbegin(); it != used.rend(); it++) {
			AndNode *old = *it;
			for (auto pointee : old->pointees()) {
				if (pointee->refs) {
					printf("%p", pointee);
					throw std::runtime_error("Loops found");
				}
			}
		}

		std::vector<AndNode> new_storage;
		new_storage.resize(used.size());
		int pos = 0;
		for (auto it = used.rbegin(); it != used.rend(); it++) {
			AndNode *old = *it;
			AndNode *new_ = &new_storage[pos++];
			old->apply_replacements();
			*new_ = *old;
			old->replacement = new_;
		}

		int nremoved = node_storage.size() - new_storage.size();
		new_storage.swap(node_storage);

		if (nremoved)
			printf("Removed %d unused nodes\n", nremoved);
		return nremoved;
	}

	int frontier()
	{
		int frontier_size = 1; // first item is special (used for PO scratch)
		std::vector<int> free_indices;

		for (auto node : nodes) {
			node->fid = 0;
			node->visited = false;
		}

		for (auto it = nodes.rbegin(); it != nodes.rend(); ++it) {
			for (auto node_repr : (*it)->fanins()) {
				for (AndNode *node = node_repr;
						node != NULL; node = node->sibling) {
					assert(!node->visited);
					if (!node->fid) {
						if (free_indices.empty())
							free_indices.push_back(frontier_size++);
						node->fid = free_indices.back();
						free_indices.pop_back();
					}
				}
			}

			// free our index
			if ((*it)->fid != 0)
				free_indices.push_back((*it)->fid);
			(*it)->visited = true;
		}

		printf("Frontier is %d wide at its peak\n", frontier_size);
		return frontier_size;
	}

	void invalidate_matches()
	{
		matches_valid = false;
		match_storage.clear();
	}

	static bool cut_union(AndNode *target[], int &cutlen, int max_cut, CutList in1, CutList in2)
	{
		CutList::iterator it2 = in2.begin();

		for (auto n1 : in1) {
			for (; it2 != in2.end() && *it2 < n1; ++it2) {
				if (cutlen == max_cut)
					return false;
				target[cutlen++] = *it2;
			}
			if (it2 != in2.end() && *it2 == n1)
				++it2;
			if (cutlen == max_cut)
				return false;
			target[cutlen++] = n1;
		}

		for (; it2 != in2.end(); ++it2) {
			if (cutlen == max_cut)
				return false;
			target[cutlen++] = *it2;
		}
		return true;
	}

	void prepare_cuts(int npriority_cuts, int nmatches_max, int max_cut=CUT_MAXIMUM)
	{
		if (max_cut < 3 || max_cut > CUT_MAXIMUM)
			throw std::runtime_error("Maximum cut size out of range");

		if (npriority_cuts < 1 || npriority_cuts > 65536)
			throw std::runtime_error("Priority cuts number out of range");

		invalidate_matches();
		int frontier_size = frontier();

		struct PriorityCut {
			AndNode *cut[CUT_MAXIMUM];
			truth6 function;
#ifdef SIBLING_RECORDING
			int nused_siblings = 0;
			AndNode *used_siblings[16];

			void insert_sibling(AndNode *sibl)
			{
				int i;
				for (i = 0; i < nused_siblings; i++)
					if (used_siblings[i] >= sibl)
						break;

				if (i < nused_siblings && used_siblings[i] == sibl)
					return;

				nused_siblings++;
				assert(nused_siblings <= std::size(used_siblings));
				for (; i < nused_siblings; i++)
					std::swap(sibl, used_siblings[i]);
			}
#endif
		};
		struct NodeCache {
			int ps_len;
			PriorityCut *ps;
			AndNode *mark;
		};
		PriorityCut *pcuts = new PriorityCut[frontier_size * npriority_cuts];
		NodeCache *cache = new NodeCache[frontier_size];

		int matches_remaining = 0;
		size_t matches_allocated = 0;
		AndNode::Match *matches_page = NULL;

		int nnodes = 0, nsatur_cuts = 0, nsatur_matches = 0;

		uint64_t nmatches_sum = 0;
		uint64_t nmatches_sum_geom = 0;

		// Go over the nodes in topological order
		for (auto node : nodes.w_indication()) {
			if (matches_remaining < (nmatches_max + 1)) {
				matches_remaining = nmatches_max * 128;
				matches_allocated += matches_remaining * sizeof(AndNode::Match);
				matches_page = new AndNode::Match[matches_remaining];
				match_storage.emplace_back(matches_page);
			}
			node->matches = matches_page;

			// Clear the cache
			NodeCache *lcache = &cache[node->fid];
			lcache->ps = &pcuts[node->fid * npriority_cuts];
			lcache->ps_len = 0;
			lcache->mark = node;

			if (node->pi)
				continue;

			if (node->po) {
				// PO has empty cache
				lcache->ps_len = 0;

				// TODO: this isn't an aiger invariant, probably
				assert(node->ins[1].eval());

				auto &match = node->matches[0];
				node->matches[1].cut[0] = NULL;
				node->pol[0].sel = 0;

				int cutlen = 0;
				for (auto fanin : node->fanins())
					match.cut[cutlen++] = fanin;
				match.cut[cutlen++] = NULL;

				if (!node->ins[0].node) {
					match.semiclass = npn_semiclass(0, 0, match.npn);
					static Target dummy = {
						.cell = NULL,
						.map = NPN::identity(0),
					};
					node->pol[0].sel_target = &dummy;
				} else {
					static Target dummy = {
						.cell = NULL,
						.map = NPN::identity(1),
					};
					node->pol[0].sel_target = &dummy;
					match.semiclass = npn_semiclass(node->ins[0].negated ? 1 : 2, 1,
													match.npn);
				}
				matches_remaining -= 2;
				matches_page += 2;
				continue;
			}

			AndNode *n1 = node->ins[0].node, *n1_save = n1;
			AndNode *n2 = node->ins[1].node;

			std::set<int> seen_cuts;

			// Find up to nmatches_max of matches to technology cells
			int nmatches = 0;

			bool n1_negated = node->ins[0].negated;
			bool n2_negated = node->ins[1].negated;

			bool n1_choicing = false, n2_choicing = false;
		choice_switched:
			assert(n1 && n2);
			assert(cache[n1->fid].mark == n1);
			assert(cache[n2->fid].mark == n2);

			AndNode *t1_nodes[2] = { n1, NULL };
			AndNode *t2_nodes[2] = { n2, NULL };
			CutList t1(t1_nodes);
			CutList t2(t2_nodes);

			for (int i = -1; i < cache[n1->fid].ps_len; i++)
			for (int j = -1; j < cache[n2->fid].ps_len; j++) {
				CutList n1_cut = ((i == -1) ? t1 : cache[n1->fid].ps[i].cut);
				CutList n2_cut = ((j == -1) ? t2 : cache[n2->fid].ps[j].cut);

				truth6 n1_function = ((i == -1) ? 2 : cache[n1->fid].ps[i].function);
				if (n1_negated) n1_function ^= mask6(n1_cut.size);
				truth6 n2_function = ((j == -1) ? 2 : cache[n2->fid].ps[j].function);
				if (n2_negated) n2_function ^= mask6(n2_cut.size);

				AndNode *working_cut[CUT_MAXIMUM];
				int cutlen = 0;
				if (!cut_union(working_cut, cutlen, max_cut, n1_cut, n2_cut))
					continue;

				truth6 cut_function = recode6(n1_function, n1_cut, working_cut) &
										recode6(n2_function, n2_cut, working_cut);
				{
					uint32_t removal_mask;
					cut_function = reduce6(cut_function, cutlen, removal_mask);
					int cutlen2 = 0;
					for (int m = 0; m < cutlen; m++) {
						if (!(removal_mask & 1 << m))
							working_cut[cutlen2++] = working_cut[m];
					}
					cutlen = cutlen2;
					if (cutlen < CUT_MAXIMUM)
						working_cut[cutlen] = NULL;
				}

				int hash = 0;
				for (auto node : CutList{working_cut})
					hash = ((hash << 5) + hash) ^ (uintptr_t) node;
				if (seen_cuts.count(hash))
					continue;
				seen_cuts.insert(hash);

				if (lcache->ps_len == npriority_cuts)
					continue;
				int slot = lcache->ps_len++;

				std::copy(working_cut, working_cut + CUT_MAXIMUM, lcache->ps[slot].cut);
				lcache->ps[slot].function = cut_function;
#ifdef SIBLING_RECORDING
				{
					auto &ps = lcache->ps[slot];
					ps.nused_siblings = 0;
					if (n1_choicing)
						ps.insert_sibling(n1);
					if (n2_choicing)
						ps.insert_sibling(n2);
					if (i >= 0) {
						auto &ps1 = cache[n1->fid].ps[i];
						for (int k = 0; k < ps1.nused_siblings; k++)
							ps.insert_sibling(ps1.used_siblings[k]);
					}
					if (j >= 0) {
						auto &ps2 = cache[n2->fid].ps[i];
						for (int k = 0; k < ps2.nused_siblings; k++)
							ps.insert_sibling(ps2.used_siblings[k]);
					}
				}
#endif

				NPN npn;
				truth6 semiclass = npn_semiclass(cut_function, cutlen, npn);
				if (target_index.classes.count(std::make_pair(semiclass, cutlen)) && nmatches < nmatches_max) {
					auto &match = node->matches[nmatches++];
					match.semiclass = semiclass;
					match.npn = npn;
					std::copy(working_cut, working_cut + CUT_MAXIMUM,
							  match.cut);
#ifdef SIBLING_RECORDING
					{
						auto &ps = lcache->ps[slot];
						std::copy(ps.used_siblings,
								  ps.used_siblings + ps.nused_siblings,
								  std::back_inserter(match.used_siblings));
					}
#endif
				}
			}
			if (n1->sibling) {
				n1_negated ^= n1->polarity ^ n1->sibling->polarity;
				n1 = n1->sibling;
				n1_choicing = true;
				goto choice_switched;
			} 
			if (n2->sibling) {
				n2_negated ^= n2->polarity ^ n2->sibling->polarity;
				n2 = n2->sibling;
				n1_negated ^= n1->polarity ^ n1_save->polarity;
				n1 = n1_save;
				n1_choicing = false;
				n2_choicing = true;
				goto choice_switched;
			} 

			node->matches[nmatches].cut[0] = NULL;
			matches_remaining -= nmatches + 1;
			matches_page += nmatches + 1;

			nnodes++;
			if (nmatches == nmatches_max)
				nsatur_matches++;
			if (lcache->ps_len == npriority_cuts)
				nsatur_cuts++;

			nmatches_sum += nmatches;
			nmatches_sum_geom += (uint64_t) nmatches * nmatches;
		}

		delete[] pcuts;
		delete[] cache;

		printf("\nCut matching statistics:\n");
		printf("  %d nodes", nnodes);
		printf(" %4.2f MiB cut cache", ((float) frontier_size * npriority_cuts) / (1024 * 1024));
		printf(" %4.2f MiB match cache\n", ((float) matches_allocated) / (1024 * 1024));
		printf("  saturated %d cuts (%.1f %%),", nsatur_cuts, ((float) nsatur_cuts * 100) / nnodes);
		printf(" %d matches (%.1f %%)\n", nsatur_matches, ((float) nsatur_matches * 100) / nnodes);
		printf("  matches %.1f mean %.1f geom\n", (float) nmatches_sum / nnodes,
			   sqrt((float) nmatches_sum_geom / nnodes));
		printf("\n");

		// no mapping on top of the matches yet
		for (auto node : nodes)
		for (int C = 0; C < 2; C++) {
			node->pol[C].map_fouts = 0;
		}

		if (target_index.inv_cell && target_index.tie.cell)
			matches_valid = true;
		else
			printf("Missing basic inverter/tie high/tie low cells\n");
	}

	void lose_choices()
	{
		// touching the siblings invalidates matches
		invalidate_matches();

		int nsiblings = 0;
		for (auto node : nodes) {
			if (node->sibling)
				nsiblings++;
			node->sibling = NULL;
		}
		printf("Cleared %d choice pairs\n", nsiblings);
	}

	static void deref_cut(AndNode *node, bool C)
	{
		if (node->pi)
			return;

		int n = 0;
		auto &match = node->matches[node->pol[C].sel];
		assert(node->pol[C].sel_target);
		NPN local_map = node->pol[C].sel_target->map * match.npn;
		for (auto cut_node : CutList{match.cut}) {
			bool cut_nodeC = local_map.ic[n++];
			assert(cut_node != node);
			auto &map_fouts = cut_node->pol[cut_nodeC].map_fouts;
			assert(map_fouts >= 1);
			if (!--map_fouts)
				deref_cut(cut_node, cut_nodeC);
		}
	}

	static float ref_cut(AndNode *node, bool C)
	{
		if (node->pi)
			return 0;

		float sum = 0;
		int n = 0;
		auto &match = node->matches[node->pol[C].sel];
		assert(node->pol[C].sel_target);
		NPN local_map = node->pol[C].sel_target->map * match.npn;
		for (auto cut_node : CutList{match.cut}) {
			assert(!cut_node->po);
			bool cut_nodeC = local_map.ic[n++];
			assert(cut_node != node);
			auto &cut_pol = cut_node->pol[cut_nodeC];
			if (!cut_pol.map_fouts++) {
				if (cut_node->pi && cut_nodeC) {
					sum += target_index.inv_cell->area();
				} else if (!cut_node->pi) {
					assert(cut_pol.sel_target);
					sum += cut_pol.sel_target->cell->area() + ref_cut(cut_node, cut_nodeC);
				}
			}
			assert(cut_pol.map_fouts >= 1);
		}
		return sum;
	}

	void ensure_matches()
	{
		if (!matches_valid)
			throw std::runtime_error("Matches are not pre-computed, run 'prepare_cuts' first");
	}

	float walk_mapping()
	{
		ensure_matches();

		for (auto it = nodes.rbegin(); it != nodes.rend(); ++it) {
			AndNode *node = *it;
			assert(node->pol[1].map_fouts == 0);
			assert(node->pol[0].map_fouts >= 0);
			assert(node->pol[0].map_fouts <= node->po ? 1 : 0);

			if (node->pol[0].map_fouts)
				deref_cut(node, false);
			node->pol[0].map_fouts = 0;
		}

		for (auto node : nodes)
			assert(!node->pol[0].map_fouts && !node->pol[1].map_fouts);

		float area = 0;

		for (auto node : nodes)
		if (node->po) {
			assert(node->pol[0].map_fouts == 0);
			assert(node->pol[1].map_fouts == 0);
			if (!node->pol[0].map_fouts++)
				area += ref_cut(node, false);
		}

		return area;
	}

	void exact_round()
	{
		ensure_matches();

		for (auto node : nodes) {
			if (node->pi || node->po)
				continue;

			for (int C = 0; C < 2; C++) {
				auto &pol = node->pol[C];

				if (pol.map_fouts)
					deref_cut(node, C);

				float best_area = std::numeric_limits<float>::max();
				int best_index = -1;
				Target *best_target = nullptr;

				for (int i = 0;; i++) {
					auto &match = node->matches[i];
					if (!match.cut[0])
						break;

					for (auto &target : target_index.classes.at(std::make_pair(match.semiclass, CutList{match.cut}.size))) {
						// make sure this is a valid match
						if ((target.map * match.npn).oc != C)
							continue;

						pol.sel = i;
						pol.sel_target = &target;
						float area = target.cell->area() + ref_cut(node, C);
						if (area < best_area) {
							best_area = area;
							best_index = i;
							best_target = &target;
						}
						deref_cut(node, C);
					}
				}

				assert(best_index != -1);
				pol.sel = best_index;
				pol.sel_target = best_target;

				if (pol.map_fouts)
					ref_cut(node, C);
			}
		}
	}

	void annealing_round(float temp)
	{
		ensure_matches();

		std::random_device rd;
    	std::mt19937 gen(rd());
		std::exponential_distribution<> d(1.0f / temp);

		for (auto node : nodes) {
			if (node->pi || node->po)
				continue;

			for (int C = 0; C < 2; C++) {
				auto &pol = node->pol[C];

				if (pol.map_fouts)
					deref_cut(node, C);

				float best_area = std::numeric_limits<float>::max();
				int best_index = -1;
				Target *best_target = nullptr;

				for (int i = 0;; i++) {
					auto &match = node->matches[i];
					if (!match.cut[0])
						break;

					for (auto &target : target_index.classes.at(std::make_pair(match.semiclass, CutList{match.cut}.size))) {
						// make sure this is a valid match
						if ((target.map * match.npn).oc != C)
							continue;

						pol.sel = i;
						pol.sel_target = &target;
						float area = target.cell->area() + ref_cut(node, C) \
										- d(gen);
						if (area < best_area) {
							best_area = area;
							best_index = i;
							best_target = &target;
						}
						deref_cut(node, C);
					}
				}

				assert(best_index != -1);
				pol.sel = best_index;
				pol.sel_target = best_target;

				if (pol.map_fouts)
					ref_cut(node, C);
			}
		}
	}

	void depth_round(bool first)
	{
		ensure_matches();
		fanouts(true);

		for (auto node : nodes) {
			if (node->po)
				continue;

			if (node->pi) {
				int fanouts = first ? node->fanouts : std::max(node->pol[1].map_fouts, 1);
				node->pol[0].farea = 0;
				node->pol[1].farea = target_index.inv_cell->area() / fanouts;
				node->pol[0].depth = 0;
				node->pol[1].depth = 1;
				continue;
			}

			for (int C = 0; C < 2; C++) {
				auto &pol = node->pol[C];

				float best_area = std::numeric_limits<float>::max();
				int best_index = -1;
				int best_depth = std::numeric_limits<int>::max();
				Target *best_target = nullptr;

				int fanouts = first ? node->fanouts : std::max(pol.map_fouts, 1);

				for (int i = 0;; i++) {
					auto &match = node->matches[i];
					if (!match.cut[0])
						break;

					for (auto &target : target_index.classes.at(std::make_pair(match.semiclass, CutList{match.cut}.size))) {
						// make sure this is a valid match
						if ((target.map * match.npn).oc != C)
							continue;

						NPN local_map = target.map * match.npn;
						float area = target.cell->area();
						int depth = 0;
						int n = 0;
						for (auto cut_node : CutList{match.cut}) {
							bool cut_nodeC = local_map.ic[n++];
							auto &cut_pol = cut_node->pol[cut_nodeC];
							area += cut_pol.farea;
							depth = std::max(depth, cut_pol.depth + 1);
						}

						if (depth < best_depth || (depth == best_depth && area < best_area)) {
							best_depth = depth;
							best_area = area;
							best_index = i;
							best_target = &target;
						}
					}
				}

				assert(best_index != -1);

				if (pol.map_fouts && (pol.sel != best_index || pol.sel_target != best_target)) {
					deref_cut(node, C);
					pol.sel = best_index;
					pol.sel_target = best_target;
					ref_cut(node, C);
				} else {
					pol.sel = best_index;
					pol.sel_target = best_target;
				}

				pol.depth = best_depth;
				pol.farea = best_area / fanouts;
			}
		}
	}

	void depth2_round(bool first)
	{
		ensure_matches();
		fanouts(true);

		for (auto node : nodes) {
			if (node->po)
				continue;

			if (node->pi) {
				int fanouts = first ? node->fanouts : std::max(node->pol[1].map_fouts, 1);
				node->pol[0].farea = 0;
				node->pol[1].farea = target_index.inv_cell->area() / fanouts;
				node->pol[0].depth = 0;
				node->pol[1].depth = 1;
				continue;
			}

			for (int C = 0; C < 2; C++) {
				auto &pol = node->pol[C];

				float best_area = std::numeric_limits<float>::max();
				int best_index = -1;
				int best_depth = std::numeric_limits<int>::max();
				Target *best_target = nullptr;

				int fanouts = first ? node->fanouts : std::max(pol.map_fouts, 1);

				for (int i = 0;; i++) {
					auto &match = node->matches[i];
					if (!match.cut[0])
						break;

					for (auto &target : target_index.classes.at(std::make_pair(match.semiclass, CutList{match.cut}.size))) {
						// make sure this is a valid match
						if ((target.map * match.npn).oc != C)
							continue;

						NPN local_map = target.map * match.npn;
						float area = target.cell->area();
						int depth = 0;
						int n = 0;
						for (auto cut_node : CutList{match.cut}) {
							bool cut_nodeC = local_map.ic[n++];
							auto &cut_pol = cut_node->pol[cut_nodeC];
							area += cut_pol.farea;
							depth += cut_pol.depth;
						}

						if (depth < best_depth || (depth == best_depth && area < best_area)) {
							best_depth = depth;
							best_area = area;
							best_index = i;
							best_target = &target;
						}
					}
				}

				assert(best_index != -1);

				if (pol.map_fouts && (pol.sel != best_index || pol.sel_target != best_target)) {
					deref_cut(node, C);
					pol.sel = best_index;
					pol.sel_target = best_target;
					ref_cut(node, C);
				} else {
					pol.sel = best_index;
					pol.sel_target = best_target;
				}

				pol.depth = best_depth;
				pol.farea = best_area / fanouts;
			}
		}
	}

	void area_flow_round(float refs_blend)
	{
		ensure_matches();
		fanouts(true, true);

		for (auto node : nodes)
		for (int C = 0; C < 2; C++) {
			node->pol[C].flow_fouts = \
				std::max(refs_blend * node->fanouts + (1.0f - refs_blend) * node->pol[C].map_fouts, 1.0f);
		}

		for (auto node : nodes) {
			if (node->po)
				continue;

			if (node->pi) {
				node->pol[0].farea = 0;
				node->pol[1].farea = target_index.inv_cell->area() / node->pol[1].flow_fouts;
				continue;
			}

			for (int C = 0; C < 2; C++) {
				auto &pol = node->pol[C];

				float best_area = std::numeric_limits<float>::max();
				int best_index = -1;
				Target *best_target = nullptr;

				for (int i = 0; true; i++) {
					auto &match = node->matches[i];
					if (!match.cut[0])
						break;

					for (auto &target : target_index.classes.at(std::make_pair(match.semiclass, CutList{match.cut}.size))) {
						// make sure this is a valid match
						if ((target.map * match.npn).oc != C)
							continue;

						NPN local_map = target.map * match.npn;
						float area = target.cell->area();
						int n = 0;
						for (auto cut_node : CutList{match.cut}) {
							bool cut_nodeC = local_map.ic[n++];
							auto &cut_pol = cut_node->pol[cut_nodeC];
							area += cut_pol.farea;
						}
						area = std::min(area, 1e32f);

						if (area < best_area) {
							best_area = area;
							best_index = i;
							best_target = &target;
						}
					}
				}

				assert(best_index != -1);

				if (pol.map_fouts && (pol.sel != best_index || pol.sel_target != best_target)) {
					deref_cut(node, C);
					pol.sel = best_index;
					pol.sel_target = best_target;
					ref_cut(node, C);
				} else {
					pol.sel = best_index;
					pol.sel_target = best_target;
				}

				pol.area = best_area;
				pol.farea = best_area / node->pol[C].flow_fouts;
			}
		}
	}

	void fuzzy_round(bool first, float temp)
	{
		ensure_matches();
		fanouts(true, true);

		for (auto node : nodes)
		for (int C = 0; C < 2; C++) {
			if (first)
				node->pol[C].flow_fouts = std::max(node->pol[C].map_fouts, 1);
			else
				node->pol[C].flow_fouts = std::max(node->pol[C].fuzzy_fouts, 1.0f);
			node->pol[C].fuzzy_fouts = 0;
		}

		for (auto node : nodes) {
			if (node->po)
				continue;

			if (node->pi) {
				node->pol[0].farea = 0;
				node->pol[1].farea = target_index.inv_cell->area() / node->pol[1].flow_fouts;
				continue;
			}

			for (int C = 0; C < 2; C++) {
				auto &pol = node->pol[C];

				float best_area = std::numeric_limits<float>::max();
				int best_index = -1;
				Target *best_target = nullptr;
				float Z = 0.0;
				pol.area = 0;

				for (int i = 0; true; i++) {
					auto &match = node->matches[i];
					if (!match.cut[0])
						break;

					for (auto &target : target_index.classes.at(std::make_pair(match.semiclass, CutList{match.cut}.size))) {
						// make sure this is a valid match
						if ((target.map * match.npn).oc != C)
							continue;

						NPN local_map = target.map * match.npn;
						float area = target.cell->area();
						int n = 0;
						for (auto cut_node : CutList{match.cut}) {
							bool cut_nodeC = local_map.ic[n++];
							auto &cut_pol = cut_node->pol[cut_nodeC];
							area += cut_pol.farea;
						}
						area = std::min(area, 1e32f);

						if (area < best_area) {
							if (best_index != -1)
								Z *= std::exp((area - best_area) / temp);
							best_area = area;
							best_index = i;
							best_target = &target;
						}

						Z += std::exp((best_area - area) / temp);
					}
				}

				for (int i = 0; true; i++) {
					auto &match = node->matches[i];
					if (!match.cut[0])
						break;

					for (auto &target : target_index.classes.at(std::make_pair(match.semiclass, CutList{match.cut}.size))) {
						// make sure this is a valid match
						if ((target.map * match.npn).oc != C)
							continue;

						NPN local_map = target.map * match.npn;
						float area = target.cell->area();
						int n = 0;
						for (auto cut_node : CutList{match.cut}) {
							bool cut_nodeC = local_map.ic[n++];
							auto &cut_pol = cut_node->pol[cut_nodeC];
							area += cut_pol.farea;
						}
						area = std::min(area, 1e32f);

						pol.area += area * (std::exp((best_area - area) / temp) / Z);
					}
				}

				assert(best_index != -1);
				if (pol.map_fouts && (pol.sel != best_index || pol.sel_target != best_target)) {
					deref_cut(node, C);
					pol.sel = best_index;
					pol.sel_target = best_target;
					ref_cut(node, C);
				} else {
					pol.sel = best_index;
					pol.sel_target = best_target;
				}

				pol.farea = pol.area / node->pol[C].flow_fouts;
			}
		}

		for (auto it = nodes.rbegin(); it != nodes.rend(); ++it) {
			AndNode *node = *it;

			if (node->po) {
				int C = 0;
				for (int i = 0; true; i++) {
					auto &match = node->matches[i];
					if (!match.cut[0])
						break;

					for (auto &target : target_index.classes.at(std::make_pair(match.semiclass, CutList{match.cut}.size))) {
						// make sure this is a valid match
						if ((target.map * match.npn).oc != C)
							continue;

						NPN local_map = target.map * match.npn;
						float area = target.cell->area();
						int n = 0;
						for (auto cut_node : CutList{match.cut}) {
							bool cut_nodeC = local_map.ic[n++];
							auto &cut_pol = cut_node->pol[cut_nodeC];
							cut_pol.fuzzy_fouts += 1.0f;
						}
					}
				}
				continue;
			}

			if (node->pi)
				continue;

			for (int C = 0; C < 2; C++) {
				auto &pol = node->pol[C];

				float best_area = std::numeric_limits<float>::max();
				int best_index = -1;
				Target *best_target = nullptr;
				float Z = 0.0;

				for (int i = 0; true; i++) {
					auto &match = node->matches[i];
					if (!match.cut[0])
						break;

					for (auto &target : target_index.classes.at(std::make_pair(match.semiclass, CutList{match.cut}.size))) {
						// make sure this is a valid match
						if ((target.map * match.npn).oc != C)
							continue;

						NPN local_map = target.map * match.npn;
						float area = target.cell->area();
						int n = 0;
						for (auto cut_node : CutList{match.cut}) {
							bool cut_nodeC = local_map.ic[n++];
							auto &cut_pol = cut_node->pol[cut_nodeC];
							area += cut_pol.farea;
						}
						area = std::min(area, 1e32f);

						if (area < best_area) {
							if (best_index != -1)
								Z *= std::exp((area - best_area) / temp);
							best_area = area;
							best_index = i;
							best_target = &target;
						}

						Z += std::exp((best_area - area) / temp);
					}
				}

				for (int i = 0; true; i++) {
					auto &match = node->matches[i];
					if (!match.cut[0])
						break;

					for (auto &target : target_index.classes.at(std::make_pair(match.semiclass, CutList{match.cut}.size))) {
						// make sure this is a valid match
						if ((target.map * match.npn).oc != C)
							continue;

						NPN local_map = target.map * match.npn;
						float area = target.cell->area();
						int n = 0;
						for (auto cut_node : CutList{match.cut}) {
							bool cut_nodeC = local_map.ic[n++];
							auto &cut_pol = cut_node->pol[cut_nodeC];
							area += cut_pol.farea;
						}
						area = std::min(area, 1e32f);

						float p = std::max(0.05f, std::min(pol.fuzzy_fouts, 0.95f)) * std::exp((best_area - area) / temp) / Z;
						assert(p >= 0.0f && p <= 1.0f);

						n = 0;
						for (auto cut_node : CutList{match.cut}) {
							bool cut_nodeC = local_map.ic[n++];
							auto &cut_pol = cut_node->pol[cut_nodeC];
							cut_pol.fuzzy_fouts += p;
						}
					}
				}
			}
		}
	}

	void save()
	{
		ensure_matches();

		// TODO: sel validity
		for (auto node : nodes)
		for (int C = 0; C < 2; C++) {
			node->pol[C].save_sel = node->pol[C].sel;
			node->pol[C].save_sel_target = node->pol[C].sel_target;
			node->pol[C].save_area = node->pol[C].area;
		}
	}

	void stitch()
	{
		ensure_matches();

		// TODO: sel validity
		for (auto node : nodes)
		for (int C = 0; C < 2; C++) {
			if (node->pi || node->po)
				continue;

			if (node->pol[C].save_area < node->pol[C].area) {
				if (node->pol[C].map_fouts)
					deref_cut(node, C);
				node->pol[C].sel = node->pol[C].save_sel;
				node->pol[C].sel_target = node->pol[C].save_sel_target;
				if (node->pol[C].map_fouts)
					ref_cut(node, C);
			}
		}
	}

	void extract_mapping(sta::ConcreteNetwork *stan)
	{
		ensure_matches();

		sta::Library *lib = stan->findLibrary("mapping");
		sta::Cell *top_cell = lib ? stan->findCell(lib, name.c_str()) : NULL;
		if (!top_cell)
			throw std::runtime_error("Cell for the mapped network missing");

		sta::Instance *top = stan->makeInstance(top_cell, "", NULL);
		sta::Net *zero = stan->makeNet("$zero", top);
		sta::Net *one = stan->makeNet("$one", top);
		sta::Instance *hilo_cell = stan->makeInstance(target_index.tie.cell, "hilo_cell", top);
		stan->connect(hilo_cell, target_index.tie.lo, zero);
		stan->connect(hilo_cell, target_index.tie.hi, one);

		// TODO: enforce ports split per-bit
		int autoidx = 0;
		char buf[128];
		for (auto node : nodes)
		if (node->pi || node->po) {
			std::string port_name = sta::portVerilogToSta(node->label.c_str());
			sta::Port *port = stan->findPort(top_cell, port_name.c_str());
			assert(!stan->findPin(top, port));
			std::string net_name = sta::netVerilogToSta(node->label.c_str());
			node->net = stan->makeNet(net_name.c_str(), top);
			sta::Pin *pin = stan->makePin(top, port, NULL);
			stan->makeTerm(pin, node->net);

			if (node->pi) {
				node->pol[0].net = node->net;
				snprintf(buf, sizeof(buf), "_%08d_", autoidx++);
				node->pol[1].net = stan->makeNet(buf, top);

				snprintf(buf, sizeof(buf), "_%08d_", autoidx++);
				sta::Instance *inv = stan->makeInstance(target_index.inv_cell, buf, top);
				sta::LibertyPort *in, *out;
				target_index.inv_cell->bufferPorts(in, out);
				assert(in && out);
				stan->connect(inv, in, node->net);
				stan->connect(inv, out, node->pol[1].net);
			}
		}

		stan->setTopInstance(top);

		for (auto node : nodes)
		for (int C = 0; C < 2; C++) {
			if (!node->pol[C].map_fouts || node->pi || node->po)
				continue;

			auto &match = node->matches[node->pol[C].sel];
			auto target = node->pol[C].sel_target;
			NPN local_map = target->map * match.npn;
			sta::LibertyCell *cell = target->cell;
			assert(cell);

			snprintf(buf, sizeof(buf), "_%08d_", autoidx);
			sta::Instance *gate = stan->makeInstance(cell, buf, top);

			std::vector<sta::LibertyPort *> inports;
			sta::LibertyPort *outport = nullptr;
			sta::LibertyCellPortIterator it(cell);
			while (it.hasNext()) {
				sta::LibertyPort *port = it.next();
				if (port->direction()->isInput()) {
					inports.push_back(port);
				} else if (port->direction()->isOutput()) {
					outport = port;
				} else {
					abort();
				}
			}
			assert(outport);
			assert(local_map.ninputs() == (int) inports.size());
			int cutidx = 0;
			for (auto cut_node : CutList{match.cut}) {
				assert(cut_node->pol[local_map.ic[cutidx]].net);
				stan->connect(gate, inports[local_map.p[cutidx]], cut_node->pol[local_map.ic[cutidx]].net);
				cutidx++;
			}

			snprintf(buf, sizeof(buf), "_%08d_%s", autoidx++, outport->name());
			node->pol[C].net = stan->makeNet(buf, top);
			stan->connect(gate, outport, node->pol[C].net);
		}

		// only visit POs with valid mapping (those will have map_fouts set)
		for (auto node : nodes)
		if (node->po && node->pol[0].map_fouts) {
			auto &match = node->matches[node->pol[0].sel];
			CutList cut_list{match.cut};
			if (!cut_list.size) {
				if (node->ins[0].eval() && node->ins[1].eval())
					stan->mergeInto(node->net, one);
				else
					stan->mergeInto(node->net, zero);
			} else {
				assert(node->ins[0].node && !node->ins[1].node);
				assert(node->ins[1].eval());
				assert(cut_list.size == 1);
				assert(match.cut[0]->pol[node->ins[0].negated].net);
				stan->mergeInto(node->net, match.cut[0]->pol[node->ins[0].negated].net);
			}
		}
	}
};

extern "C" {
extern int Pressmold_swig_Init(Tcl_Interp *interp);
}

namespace sta {
extern const char* tcl_inits[];
extern const char *pressmold_swig_tcl_inits[];
}

#include "commands.h"

Network net;

void prepare_cuts_cmd(int cuts, int matches, int max_cut)
{
	if (max_cut == -1)
		max_cut = CUT_MAXIMUM;

	net.prepare_cuts(cuts, matches, max_cut);
}

void mapping_round_cmd(const char *kind, float param, bool param2)
{
	if (!strcmp(kind, "exact")) {
		net.exact_round();
	} else if (!strcmp(kind, "depth")) {
		net.depth_round(param != 0);
	} else if (!strcmp(kind, "depth2")) {
		net.depth2_round(param != 0);
	} else if (!strcmp(kind, "flow")) {
		net.area_flow_round(param);
	} else if (!strcmp(kind, "anneal")) {
		net.annealing_round(param);
	} else if (!strcmp(kind, "fuzzy")) {
		net.fuzzy_round(param2, param);
	} else if (!strcmp(kind, "save")) {
		net.save();
	} else if (!strcmp(kind, "stitch")) {
		net.stitch();
	} else {
		throw std::runtime_error("Unknown mapping round kind");
	}

	int nedges = 0;
	int ncells = 0;
	float area = net.walk_mapping();

	for (auto node : net.nodes)
	for (int C = 0; C < 2; C++) {
		if (node->pi || node->po || !node->pol[C].map_fouts)
			continue;
		ncells++;
		nedges += CutList{node->matches[node->pol[C].sel].cut}.size;
	}

	printf("%6s  A=%8.1f  N=%5d  E=%5d", kind, area, ncells, nedges);
	if (!strcmp(kind, "flow")) {
		printf("  (blend=%1.3f)", param);
	} else if (!strcmp(kind, "anneal")) {
		printf("  (T=%1.3f)", param);
	} else if (!strcmp(kind, "fuzzy")) {
		printf("  (T=%1.3f)%s", (float) param, param2 ? " init" : "");
	}
	printf("\n");
}

void report_mapping()
{
	net.ensure_matches();

	std::map<sta::LibertyCell *, int> cell_number;
	for (auto node : net.nodes)
	for (int C = 0; C < 2; C++)
	if (!node->pi && !node->po && node->pol[C].map_fouts) {
		assert(node->pol[C].sel_target);
		cell_number[node->pol[C].sel_target->cell]++;
	}

	printf("\n");
	int cell_no_sum = 0;
	double area_sum = 0;
	for (auto pair : cell_number) {
		cell_no_sum += pair.second;
		area_sum += pair.first->area() * pair.second;
		printf("  %26s  %4d  %.3e\n", pair.first->name(), pair.second, pair.first->area() * pair.second);
	}
	printf("\nSum: %d cells %.3e area\n\n", cell_no_sum, area_sum);
}

void report_sibling_usage()
{	
#ifdef SIBLING_RECORDING
	printf("\n");

	net.fanouts(false, true);

	for (auto node : net.nodes)
		node->visited = false;

	for (auto node : net.nodes)
	for (int C = 0; C < 2; C++)
	if (!node->pi && node->pol[C].map_fouts) {
		assert(node->pol[C].sel != -1);
		for (auto sibling : node->matches[node->pol[C].sel].used_siblings) {
			assert(!sibling->in_repr);
			sibling->visited = true;
		}
	}

	for (auto node : net.nodes)
	if (node->in_repr) {
		for (AndNode *sibl = node->sibling; sibl; sibl = sibl->sibling)
		if (sibl->visited) {
			std::set<AndNode *> support;
			std::set<AndNode *> seen = {sibl};
			std::vector<AndNode *> queue = {sibl};

			while (!queue.empty()) {
				AndNode *n = queue.back(); queue.pop_back();
				for (auto fanin : n->fanins()) {
					if (fanin->in_repr) {
						support.insert(fanin);
					} else if (!seen.count(fanin)) {
						seen.insert(fanin);
						queue.push_back(fanin);
					}
				}
			}

			std::vector<AndNode *> supp;
			std::copy(support.begin(), support.end(),
					  std::back_inserter(supp));

			auto t1 = node->truth_table(supp);
			auto t2 = sibl->truth_table(supp);
			if (sibl->polarity ^ node->polarity) {
				for (int i = 0; i < t2.size(); i++)
					t2[i] = !t2[i];
			}

			assert(!t2.empty());
			if (t1.empty()) {
				printf("support %d: spilled\n", (int) support.size());
			} else {
				int d = 0;
				for (int i = 0; i < t1.size(); i++) {
					if (t1[i] != t2[i])
						d++;
				}
				printf("support %d: diff %d\n\t", (int) support.size(), d);
				for (int i = 0; i < t1.size(); i++)
					printf("%s", t1[i] != t2[i] ? "x" : "-");
				printf("\n");
			}

			uint64_t w = 0;
			if (support.size() <= 6) {
				for (int i = 0; i < t2.size(); i++)
					if (t2[i]) w |= ((uint64_t) 1) << i;
				NPN npn;
				printf("\tsemiclass %llx (%llx)\n", npn_semiclass(w, (int) support.size(), npn), w);
			}
		}
	}
#endif
}

void extract_mapping()
{
	sta::ConcreteNetwork *stan = (sta::ConcreteNetwork *) sta::Sta::sta()->networkReader();
	net.extract_mapping(stan);
	sta::Sta::sta()->networkChanged();
}

uint64_t fexpr_eval(sta::FuncExpr &fexpr, std::vector<sta::LibertyPort *> ins)
{
	uint64_t mask = ((uint64_t) 2 << ((1 << ins.size()) - 1)) - 1;

	switch (fexpr.op()) {
	case sta::FuncExpr::op_port: {
		int i;
		for (i = 0; i < (int) ins.size(); i++)
			if (ins[i] == fexpr.port())
				break;
		assert(i != (int) ins.size());
		uint64_t ret = 0;
		for (int j = 0; j < 1 << ins.size(); j++)
			if (j & 1 << i)
				ret |= (truth6) 1 << j;
		return ret;
	}
	case sta::FuncExpr::op_not:  return mask & ~fexpr_eval(*fexpr.left(), ins);
	case sta::FuncExpr::op_or:   return fexpr_eval(*fexpr.left(), ins) | fexpr_eval(*fexpr.right(), ins);
	case sta::FuncExpr::op_and:  return fexpr_eval(*fexpr.left(), ins) & fexpr_eval(*fexpr.right(), ins);
	case sta::FuncExpr::op_xor:  return fexpr_eval(*fexpr.left(), ins) ^ fexpr_eval(*fexpr.right(), ins);
	case sta::FuncExpr::op_one:  return mask;
	case sta::FuncExpr::op_zero: return 0;
	default: abort();
	}
}

float buffer_delay(sta::Sta *sta, sta::LibertyCell *cell)
{
	sta::DcalcAnalysisPt *dcalc_ap = sta->cmdCorner()->findDcalcAnalysisPt(sta::MinMax::max());
	assert(dcalc_ap);
	sta::ArcDelayCalc *arc_delay_calc = sta->arcDelayCalc();
	assert(arc_delay_calc);
	const sta::Pvt *pvt = dcalc_ap->operatingConditions();

	sta::LibertyPort *in, *out;
	cell->bufferPorts(in, out);
	assert(in && out);

	sta::ArcDelay delay = 0;

	for (sta::TimingArcSet *arc_set : cell->timingArcSets()) {
		if (!arc_set->role()->isTimingCheck() &&
			arc_set->from() == in &&
			arc_set->to() == out)
		for (sta::TimingArc *arc : arc_set->arcs()) {
			sta::ArcDelay arc_delay;
			sta::Slew arc_slew;

			// TODO: this is wrong
			float load_cap = in->capacitance(sta::MinMax::max());

			sta::GateTimingModel *model = dynamic_cast<sta::GateTimingModel*>(arc->model());
			assert(model);
			model->gateDelay(pvt, 0.0, load_cap, false, arc_delay, arc_slew);
			model->gateDelay(pvt, arc_slew, load_cap, false, arc_delay, arc_slew);

			if (arc_delay > delay)
				delay = arc_delay;
		}
	}

	return delay;
}

bool register_cell_cmd(sta::LibertyCell *cell, bool verbose)
{
	net.matches_valid = false;

	if (cell->dontUse()) {
		if (verbose)
			printf("Ignoring dont-use cell %s\n", cell->name());
		return false;
	}

	std::vector<sta::LibertyPort *> outputs;
	std::vector<sta::LibertyPort *> inputs;

	sta::LibertyCellPortIterator it(cell);
	while (it.hasNext()) {
		sta::LibertyPort *port = it.next();
		if (port->direction()->isInput()) {
			inputs.push_back(port);
		} else if (port->direction()->isOutput()) {
			outputs.push_back(port);
		} else {
			if (verbose)
				printf("Ignoring cell %s\n", cell->name());
			return false;
		}
	}

	if (outputs.size() == 2 && inputs.size() == 0 &&
			outputs[0]->function() &&
			outputs[1]->function()) {
		sta::LibertyPort *hi = nullptr, *lo = nullptr;
		if (fexpr_eval(*outputs[0]->function(), inputs))
			hi = outputs[0];
		else
			lo = outputs[0];
		if (fexpr_eval(*outputs[1]->function(), inputs))
			hi = outputs[1];
		else
			lo = outputs[1];
		if (hi && lo) {
			auto &tie = target_index.tie;
			tie.hi = hi;
			tie.lo = lo;
			tie.cell = cell;
			if (verbose)
				printf("Detected tie cell %s\n", cell->name());
			return true;
		}
	}

	if (outputs.size() != 1 || !outputs[0]->function() || inputs.size() > 6) {
		if (verbose)
			printf("Ignoring cell %s\n", cell->name());
		return false;
	}

	truth6 print = fexpr_eval(*outputs[0]->function(), inputs);
	std::set<truth6> seen;
	if (verbose) {
		printf("Registering %s: fingerprint %llx\n inputs: ", cell->name(), print);
		for (auto port : inputs)
			printf("%s ", port->busName());
		printf("\n");
	}

	sta::Sta *sta = sta::Sta::sta();

	if (print == 0b01 && inputs.size() == 1) {
		if (!target_index.inv_cell || buffer_delay(sta, cell) < buffer_delay(sta, target_index.inv_cell))
			target_index.inv_cell = cell;
	}

	npn_semiclass_allrepr(print, inputs.size(), [&](truth6 repr, NPN &npn) {
		target_index.classes[std::make_pair(repr, inputs.size())].push_back(Target{ cell, npn.inv() });
	});

	return true;
}

// TODO: error handling
void read_aiger_cmd(const char *filename, const char *name)
{
	std::ifstream f(filename, std::ios::binary);
	if (!f.is_open())
		throw std::runtime_error(std::string("Failed to open ") + filename + "\n");
	sta::ConcreteNetwork *stan = (sta::ConcreteNetwork *) sta::Sta::sta()->networkReader();
	net = Network::read_aiger(f, stan, name);
}

void write_aig_verilog(const char *filename, const char *module_name)
{
	std::ofstream f(filename);
	if (!f.is_open())
		throw std::runtime_error(std::string("Failed to open ") + filename + "\n");
	net.dump_verilog(module_name, f);
}

void portlist_cmd()
{
	for (auto node : net.nodes) {
		if (node->pi)
			printf("input %s\n", node->label.c_str());
		if (node->po)
			printf("output %s\n", node->label.c_str());
	}
}

void report_aig()
{
	net.stats();
}

void lose_choices()
{
	net.lose_choices();
	net.consolidate();
}

static int tcl_main(Tcl_Interp *interp)
{
	int ret;

#define CHECK(f) \
	if ((ret = (f)) != TCL_OK) \
		return ret;

	CHECK(Tcl_Init(interp))
	CHECK(Tclreadline_Init(interp))
	Tcl_StaticPackage(interp, "tclreadline", Tclreadline_Init, Tclreadline_SafeInit);
	CHECK(Tcl_EvalFile(interp, TCLRL_LIBRARY "/tclreadlineInit.tcl"))

	sta::initSta();
	sta::Sta *sta = new sta::Sta;
	sta::Sta::setSta(sta);
	sta->makeComponents();
	sta->setTclInterp(interp);
	Pressmold_swig_Init(interp);
	sta::evalTclInit(interp, sta::tcl_inits);

	CHECK(Tcl_Eval(interp, "namespace import sta::*"))
	// `history` can fail. OpenROAD codebase says this is due to
	// some distributions deploying a broken script. In any case
	// purposefully don't check the error code here.
	Tcl_Eval(interp, "history");
	CHECK(Tcl_Eval(interp, "history event"))
	CHECK(Tcl_Eval(interp, "::tclreadline::readline builtincompleter true"))
	CHECK(Tcl_Eval(interp, "::tclreadline::readline customcompleter ::tclreadline::ScriptCompleter"))

	sta::evalTclInit(interp, sta::pressmold_swig_tcl_inits);
	printf(" Welcome to the Press Mold Mapper.\n Copyright 2024 Martin Povišer\n Happy pressing!\n\n");

	if (!Tcl_GetStartupScript(NULL))
		return Tcl_Eval(interp, "::tclreadline::Loop");

	return TCL_OK;
#undef CHECK
}

int main(int argc, char *argv[])
{
	Tcl_Main(argc, argv, tcl_main);
	return 0;
}
