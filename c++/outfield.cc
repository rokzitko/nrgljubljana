// outfield.cc - Code for flexible formatted output
// Copyright (C) 2009-2019 Rok Zitko

#ifndef _outfield_cc_
#define _outfield_cc_

// Storage of various statistical quantities calculated during the iteration.
class outfield;

using outfieldPtr = outfield *; // TODO: shared ptr
using vecoutptr = std::vector<outfieldPtr>;
vecoutptr allfields; // Container for all fields

// Base class for output field containers.
class outfield {
  protected:
  string _desc;     // description of the field
  string _value;    // value of the field
  double _rawvalue{}; // unformatted value stored in the field
  public:
  static int width; // width of the output field
  static int prec;  // numerical precision of the output field
  outfield() : _desc{""} {};
  // Using "pos", we can define at which position the element is inserted
  // to the list of fields. The default is at the end of the list.
  void set(const string &desc, int pos = -1) {
    _desc           = desc;
    outfieldPtr ptr = this;
    if (pos == -1)
      allfields.push_back(ptr);
    else {
      my_assert(pos < allfields.size());
      allfields.insert(begin(allfields) + pos, ptr);
    }
  };
  explicit outfield(const string &desc, int pos = -1) { set(desc, pos); }
  // We store a double precision number which is converted to an
  // appropriately formatted string.
  void setvalue(double newvalue) {
    my_assert(width != 0 && prec != 0);
    ostringstream tmp;
    tmp << setw(width) << setprecision(prec) << newvalue;
    _value    = tmp.str();
    _rawvalue = newvalue;
  }
  void operator=(double newvalue) { setvalue(newvalue); }
  string name() const { return _desc; }
  double rawvalue() const { return _rawvalue; }
  // output the header
  // Warning: *_width must be set to the correct value before
  // this member function is called!
  void putheader(ostream &F) const {
    my_assert(width != 0);
    F << setw(width) << _desc << " ";
  }
  // Output the value of the field.
  void put(ostream &F) const {
    my_assert(width != 0);
    F << setw(width) << _value << " ";
  }
};

int outfield::prec  = 0;
int outfield::width = 0;

// Return true if a field named 'query' exists.
bool outfield_exists(const string &query) {
  return any_of(begin(allfields), end(allfields), [query](const auto &i){ return i->name() == query; });
}

// Setup output fields that will appear in the file "td".
// Additional elements are defined in symmetry.cc.
namespace TD {
  outfield T, E, E2, C, F, S;
  void init() {
    allfields.clear();
    T.set("T");
    // insert others here (starting with pos=1, i.e. after the 1st element)
    E.set("<E>");
    E2.set("<E^2>");
    C.set("C");
    F.set("F");
    S.set("S");
  }
  void save_TD_quantities(ostream &F) {
    F << ' ';
    for (const auto &i : allfields) i->put(F);
    F << endl;
  }
  void save_header(ostream &F) {
    F << '#';
    for (const auto &i : allfields) i->putheader(F);
    F << endl;
  }
} // namespace TD

#endif // _outfield_cc_
