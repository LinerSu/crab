// Soundness of the tag environment used by the region domain's taint
// analysis, checked exhaustively against its concrete semantics.
//
// Concrete semantics: a state maps each location to the set of tags
// attached to its value; the empty set is "clean". Abstract element:
// separate_discrete_domain<variable, tag>, a finite map from variables
// to tag sets ordered by inclusion, where
//
//   - a key with the empty set means "clean",
//   - a key with the set T means "maybe tainted, only by tags in T",
//   - an ABSENT key means "unknown": any tag (the top of the
//     per-variable lattice). This is the only representable meaning of
//     "not mentioned", and it is what makes forget/project/top honest.
//
// gamma(A) = { sigma | forall x. sigma(x) subseteq A(x) }, gamma(bottom) = {}.
//
// The universe is three variables x, y, z; explicit tag sets range
// over the tags 1, 2, while the concrete tag universe also has a tag
// 3 that no explicit set contains (in the real domain the universe of
// tag ids is unbounded, so a finite set is never top): 512 concrete
// states, 125 abstract maps (each variable in {{}, {1}, {2}, {1,2},
// unknown}) plus bottom, with an injective gamma. Every property below
// is checked over all elements / all pairs.
#define BOOST_TEST_MODULE tag_env_soundness
#include <boost/test/unit_test.hpp>

#include "../crab_lang.hpp"

#include <crab/domains/discrete_domains.hpp>
#include <crab/domains/region/tags.hpp>
#include <crab/domains/separate_domains.hpp>
#include <crab/support/os.hpp>

#include <array>
#include <bitset>
#include <string>
#include <vector>

using namespace crab::cfg_impl;
using tag_t = crab::domains::region_domain_impl::tag<ikos::z_number>;
using env_t = crab::domains::separate_discrete_domain<z_var, tag_t>;
using set_t = env_t::mapped_type;

namespace {

constexpr unsigned NV = 3;                    // variables x, y, z
constexpr unsigned NT = 2;                    // tags in explicit sets: 1, 2
constexpr unsigned NTC = 3;                   // concrete tags: 1, 2, 3
constexpr unsigned NCONC = 1u << (NV * NTC);  // concrete states
constexpr unsigned UNKNOWN = 1u << NT;        // per-variable code for "absent"
constexpr unsigned NVAL = UNKNOWN + 1;        // per-variable codes
using gamma_t = std::bitset<NCONC>;

struct universe {
  variable_factory_t vfac;
  std::vector<z_var> vars;
  std::vector<tag_t> tags;

  universe() {
    vars.push_back(z_var(vfac["x"], crab::INT_TYPE, 32));
    vars.push_back(z_var(vfac["y"], crab::INT_TYPE, 32));
    vars.push_back(z_var(vfac["z"], crab::INT_TYPE, 32));
    tags.push_back(tag_t(ikos::z_number(1)));
    tags.push_back(tag_t(ikos::z_number(2)));
  }

  set_t mk_set(unsigned mask) const {
    set_t s = set_t::bottom();
    for (unsigned t = 0; t < NT; ++t) {
      if ((mask >> t) & 1u) {
        s += tags[t];
      }
    }
    return s;
  }

  // code < UNKNOWN: explicit tag set; code == UNKNOWN: key absent.
  env_t mk_env(const std::array<unsigned, NV> &code) const {
    env_t e;
    for (unsigned i = 0; i < NV; ++i) {
      if (code[i] < UNKNOWN) {
        e.set(vars[i], mk_set(code[i]));
      }
    }
    return e;
  }

  // The set of concrete tags a value abstracted by s may carry, as a
  // mask over the NTC concrete tags.
  unsigned mask_of(const set_t &s) const {
    if (s.is_top()) {
      return (1u << NTC) - 1;
    }
    unsigned m = 0;
    for (unsigned t = 0; t < NT; ++t) {
      if (set_t(tags[t]) <= s) {
        m |= 1u << t;
      }
    }
    return m;
  }

  gamma_t gamma(const env_t &e) const {
    gamma_t g;
    if (e.is_bottom()) {
      return g;
    }
    std::array<unsigned, NV> allowed;
    for (unsigned i = 0; i < NV; ++i) {
      allowed[i] = mask_of(e.at(vars[i]));
    }
    for (unsigned sigma = 0; sigma < NCONC; ++sigma) {
      bool ok = true;
      for (unsigned i = 0; i < NV && ok; ++i) {
        unsigned actual = (sigma >> (i * NTC)) & ((1u << NTC) - 1);
        ok = (actual & ~allowed[i]) == 0;
      }
      g[sigma] = ok;
    }
    return g;
  }

  std::vector<env_t> all() const {
    std::vector<env_t> res;
    std::array<unsigned, NV> code;
    for (code[0] = 0; code[0] < NVAL; ++code[0]) {
      for (code[1] = 0; code[1] < NVAL; ++code[1]) {
        for (code[2] = 0; code[2] < NVAL; ++code[2]) {
          res.push_back(mk_env(code));
        }
      }
    }
    res.push_back(env_t::bottom());
    return res;
  }
};

bool subset(const gamma_t &a, const gamma_t &b) { return (a & ~b).none(); }

bool same(const env_t &a, const env_t &b) { return a <= b && b <= a; }

std::string str(const env_t &e) {
  crab::crab_string_os os;
  e.write(os);
  return os.str();
}

std::string str(const set_t &s) {
  crab::crab_string_os os;
  s.write(os);
  return os.str();
}

std::string str(const z_var &v) {
  crab::crab_string_os os;
  os << v;
  return os.str();
}

} // namespace

BOOST_AUTO_TEST_CASE(order_is_gamma_inclusion) {
  universe U;
  auto elems = U.all();
  unsigned checked = 0;
  for (const env_t &a : elems) {
    gamma_t ga = U.gamma(a);
    for (const env_t &b : elems) {
      gamma_t gb = U.gamma(b);
      bool leq = a <= b;
      bool incl = subset(ga, gb);
      BOOST_CHECK_MESSAGE(leq == incl, "A=" << str(a) << " B=" << str(b)
                                            << " A<=B=" << leq
                                            << " gamma(A)<=gamma(B)=" << incl);
      ++checked;
    }
  }
  BOOST_CHECK_EQUAL(checked, elems.size() * elems.size());
}

BOOST_AUTO_TEST_CASE(join_is_sound_and_least) {
  universe U;
  auto elems = U.all();
  for (const env_t &a : elems) {
    gamma_t ga = U.gamma(a);
    for (const env_t &b : elems) {
      gamma_t gb = U.gamma(b);
      env_t j = a | b;
      gamma_t gj = U.gamma(j);
      BOOST_CHECK_MESSAGE(subset(ga | gb, gj),
                          "join unsound: " << str(a) << " | " << str(b)
                                           << " = " << str(j));
      BOOST_CHECK_MESSAGE(a <= j && b <= j,
                          "join not an upper bound: " << str(a) << " | "
                                                      << str(b) << " = "
                                                      << str(j));
      BOOST_CHECK_MESSAGE(same(j, b | a), "join not commutative: " << str(a)
                                                                    << ", "
                                                                    << str(b));
      for (const env_t &c : elems) {
        if (a <= c && b <= c) {
          BOOST_CHECK_MESSAGE(j <= c, "join not least: " << str(a) << " | "
                                                        << str(b) << " = "
                                                        << str(j)
                                                        << " but upper bound "
                                                        << str(c));
        }
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(meet_is_exact_and_greatest) {
  universe U;
  auto elems = U.all();
  for (const env_t &a : elems) {
    gamma_t ga = U.gamma(a);
    for (const env_t &b : elems) {
      gamma_t gb = U.gamma(b);
      env_t m = a & b;
      gamma_t gm = U.gamma(m);
      // sound: contains the intersection; exact: nothing more
      BOOST_CHECK_MESSAGE(subset(ga & gb, gm), "meet unsound: " << str(a)
                                                                << " & "
                                                                << str(b)
                                                                << " = "
                                                                << str(m));
      BOOST_CHECK_MESSAGE(subset(gm, ga & gb), "meet not exact: " << str(a)
                                                                  << " & "
                                                                  << str(b)
                                                                  << " = "
                                                                  << str(m));
      BOOST_CHECK_MESSAGE(m <= a && m <= b,
                          "meet not a lower bound: " << str(a) << " & "
                                                     << str(b) << " = "
                                                     << str(m));
      BOOST_CHECK_MESSAGE(same(m, b & a), "meet not commutative: " << str(a)
                                                                    << ", "
                                                                    << str(b));
      for (const env_t &c : elems) {
        if (c <= a && c <= b) {
          BOOST_CHECK_MESSAGE(c <= m, "meet not greatest: " << str(a) << " & "
                                                           << str(b) << " = "
                                                           << str(m)
                                                           << " but lower bound "
                                                           << str(c));
        }
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(top_and_bottom) {
  universe U;
  auto elems = U.all();
  env_t top;
  env_t bot = env_t::bottom();
  BOOST_CHECK(top.is_top());
  BOOST_CHECK(!top.is_bottom());
  BOOST_CHECK(bot.is_bottom());
  BOOST_CHECK(!bot.is_top());
  BOOST_CHECK(U.gamma(top).all());
  BOOST_CHECK(U.gamma(bot).none());
  for (const z_var &v : U.vars) {
    BOOST_CHECK(top.at(v).is_top());
    BOOST_CHECK(top.find(v) == nullptr);
  }
  for (const env_t &a : elems) {
    BOOST_CHECK(bot <= a);
    BOOST_CHECK(a <= top);
    BOOST_CHECK(same(a | bot, a));
    BOOST_CHECK(same(a & top, a));
    BOOST_CHECK((a & bot).is_bottom());
    BOOST_CHECK((a | top).is_top());
  }
}

// forget (operator-=) and project are existential quantification:
// the result contains every state of the operand (they move up), the
// removed variables become unknown, the others are untouched.
BOOST_AUTO_TEST_CASE(forget_and_project_move_up) {
  universe U;
  auto elems = U.all();
  for (const env_t &a : elems) {
    if (a.is_bottom()) {
      continue;
    }
    gamma_t ga = U.gamma(a);
    for (unsigned k = 0; k < NV; ++k) {
      env_t b(a);
      b -= U.vars[k];
      BOOST_CHECK_MESSAGE(subset(ga, U.gamma(b)),
                          "forget moved down: " << str(a) << " -= "
                                                << str(U.vars[k]) << " = "
                                                << str(b));
      BOOST_CHECK(a <= b);
      BOOST_CHECK(b.at(U.vars[k]).is_top());
      for (unsigned i = 0; i < NV; ++i) {
        if (i != k) {
          BOOST_CHECK(b.at(U.vars[i]) == a.at(U.vars[i]));
        }
      }
    }
    env_t p(a);
    p.project({U.vars[0]});
    BOOST_CHECK_MESSAGE(subset(ga, U.gamma(p)),
                        "project moved down: " << str(a) << " -> " << str(p));
    BOOST_CHECK(a <= p);
    BOOST_CHECK(p.at(U.vars[0]) == a.at(U.vars[0]));
    BOOST_CHECK(p.at(U.vars[1]).is_top());
    BOOST_CHECK(p.at(U.vars[2]).is_top());
  }
}

// How clean / maybe tainted / unknown are stored and read back.
BOOST_AUTO_TEST_CASE(clean_is_explicit_unknown_is_absent) {
  universe U;
  const z_var &x = U.vars[0];
  env_t e;
  // clean: an explicit empty set, kept in the map
  e.set(x, set_t::bottom());
  BOOST_CHECK(e.find(x) != nullptr);
  BOOST_CHECK(e.at(x).is_bottom());
  BOOST_CHECK(!e.is_top());
  // maybe tainted by tag 1
  e.set(x, U.mk_set(1));
  BOOST_CHECK(e.at(x) == U.mk_set(1));
  // unknown: setting top removes the key; forgetting does the same
  e.set(x, set_t::top());
  BOOST_CHECK(e.find(x) == nullptr);
  BOOST_CHECK(e.at(x).is_top());
  BOOST_CHECK(e.is_top());
  e.set(x, U.mk_set(1));
  e -= x;
  BOOST_CHECK(e.find(x) == nullptr);
  BOOST_CHECK(e.at(x).is_top());
  // rename keeps the value of an explicitly clean variable
  env_t r;
  r.set(x, set_t::bottom());
  r.rename({x}, {U.vars[1]});
  BOOST_CHECK(r.at(x).is_top());
  BOOST_CHECK(r.at(U.vars[1]).is_bottom());
  // printing: clean variables are listed, unknown ones are not
  BOOST_CHECK_EQUAL(str(env_t()), "{}");
  BOOST_CHECK_EQUAL(str(env_t::bottom()), "_|_");
  env_t w;
  w.set(x, set_t::bottom());
  w.set(U.vars[1], U.mk_set(1));
  std::string sw = str(w);
  BOOST_CHECK_MESSAGE(sw.find("x -> " + str(set_t::bottom())) != std::string::npos, sw);
  BOOST_CHECK_MESSAGE(sw.find("y -> " + str(U.mk_set(1))) != std::string::npos, sw);
  BOOST_CHECK_MESSAGE(sw.find("z") == std::string::npos, sw);
}

// The witness of docs/taint_dfa/issue.md at the lattice level: a
// variable tainted on one side of a join and not mentioned on the
// other is unknown after the join (a check on it must fail), while a
// variable explicitly clean on the other side keeps its taint.
BOOST_AUTO_TEST_CASE(one_sided_join_witness) {
  universe U;
  const z_var &x = U.vars[0];
  env_t a;
  a.set(x, U.mk_set(1));
  env_t absent;
  env_t clean;
  clean.set(x, set_t::bottom());

  BOOST_CHECK((a | absent).at(x).is_top());
  BOOST_CHECK((a | clean).at(x) == U.mk_set(1));
  BOOST_CHECK((clean | absent).at(x).is_top());
  BOOST_CHECK((a | a).at(x) == U.mk_set(1));
  env_t a2;
  a2.set(x, U.mk_set(2));
  BOOST_CHECK((a | a2).at(x) == U.mk_set(3));

  BOOST_CHECK(a <= absent);
  BOOST_CHECK(!(absent <= a));
  BOOST_CHECK(clean <= a);
  BOOST_CHECK(!(a <= clean));

  // meet: unknown adds no fact, clean wins
  BOOST_CHECK((a & absent).at(x) == U.mk_set(1));
  BOOST_CHECK((a & clean).at(x).is_bottom());
  BOOST_CHECK((a & a2).at(x).is_bottom());
}
