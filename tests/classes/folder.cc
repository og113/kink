/*
 	definitions for the Folder class and dependencies
 */

#include <string>
#include <iostream>
#include <vector>
#include <utility> // for pair
#include "simple.h"
#include "folder.h"

/*
	declarations for the FilenameAttributes class, a base class to be inherited from.
		- copy
		- copy constructor
		- operator=
		- <<
*/

// copy
void FilenameAttributes::copy(const FilenameAttributes& fa) {
	Directory 	= fa.Directory;
	Timenumber 	= fa.Timenumber;
	ID			= fa.ID;
	Suffix		= fa.Suffix;
	Extras		= fa.Extras;
}

// copy constructor
FilenameAttributes::FilenameAttributes(const FilenameAttributes& fa) {
	copy(fa);
}

//operator=
FilenameAttributes& FilenameAttributes::operator=(const FilenameAttributes& rhs) {
	copy(rhs);
	return *this;
}

ostream& operator<<(ostream& os, const FilenameAttributes& fa) {
	os << "Directory:  " << fa.Directory << endl;
	os << "Timenumber: " << fa.Timenumber << endl;
	os << "ID:         " << fa.ID << endl;
	os << "Extras:     ";
	if ((fa.Extras).size()>0) {
		for (unsigned int l=0; l<(fa.Extras).size(); l++) {
			os << (fa.Extras[l]).first << ", " << (fa.Extras[l]).second << endl;
		}
	}
	os << endl;
	os << "Suffix:     " << fa.Suffix << endl;
	return os;
}

/*
	declarations for the Filename class, publicly inherited from FilenameAttributes.
		- operator=
		- constructor(const string& filename)
		- operator()
		- <<
*/

// operator=
Filename& Filename::operator=(const Filename& rhs) {
	FilenameAttributes::operator=(rhs);
}

// constructor(const string& filename)
Filename::Filename(const string& f): FilenameAttributes() {
	string temp = f;
	size_t stop;
	stop = f.find("/");
	if (stop!=string::npos) {
		Directory = f.substr(0,stop);
		temp = temp.substr(stop+1);
	}
	if (temp.find_first_of("0123456789")==0) {
		stop = temp.find_first_not_of("0123456789");
		Timenumber = temp.substr(0,stop);
		temp = temp.substr(stop);
	}
	if (temp.find_first_not_of("_")==0) {
	 stop = temp.find_first_of("_.");
	 ID = temp.substr(0,stop);
	 temp = temp.substr(stop);
	}
	if (temp[0]=='_') {
		temp = temp.substr(1);
		while (stop!=string::npos && temp[stop]!='.') {
			stop = temp.find("_");
			if (stop==string::npos) {
				cerr << "Filename error: Extras not in pairs." << endl;
				break;
			}
			StringPair sp;
			sp.first = temp.substr(0,stop);
			temp = temp.substr(stop+1);
			stop = temp.find_first_of("_.");
			sp.second = temp.substr(0,stop);
			Extras.push_back(sp);
		}
	}
	if (stop!=string::npos && temp[stop]=='.') {
		Suffix = temp.substr(stop);
	}
}

// operator()
string Filename::operator()() const {
	string filename = Directory + "/" + Timenumber + ID;
	for (unsigned int l=0; l<Extras.size(); l++) {
		filename += "_" + Extras[l].first + "_" + Extras[l].second;
	}
	filename += Suffix;
	return filename;
}

// <<
ostream& operator<<(ostream& os, const Filename& f) {
	os << f();
	return os;
}

/*
	defintions for the FilenameComparator class, which is used by the Folder class. FilenameComparator sees the ugly details of Filename.
		- copy
		- copy constructor
		- check
		- constructor(lower,upper)
		- operator=
		- set
		- setLower
		- setUpper
		- operator(Filename)
		- <<
*/

// copy
void FilenameComparator::copy(const FilenameComparator& fc) {
	Lower = fc.Lower;
	Upper = fc.Upper;
}

// copy constructor
FilenameComparator::FilenameComparator(const FilenameComparator& fc) {
	copy(fc);
}

// check
bool FilenameComparator::check(const FilenameAttributes& l, const FilenameAttributes& u) {
	if ((l.Directory).compare(u.Directory)!=0) {
		cerr << "FilenameComparator error: Lower.Directory = " << l.Directory << ", Upper.Directory = " << u.Directory << endl;
		return false;
		}
	if ((l.ID).compare(u.ID)!=0) {
		cerr << "FilenameComparator error: Lower.ID = " << l.ID << ", Upper.ID = " << u.ID << endl;
		return false;
	}
	if ((l.Suffix).compare(u.Suffix)!=0) {
		cerr << "FilenameComparator error: Lower.Suffix = " << l.Suffix << ", Upper.Suffix = " << u.Suffix << endl;
		return false;
	if ((l.Extras).size()!=(u.Extras).size()) {
		cerr << "FilenameComparator error: Lower.Extras.size() = " << (l.Extras).size() << ", Upper.Extras.size = "\
			 << (u.Extras).size() << endl;
		return false;
	}
	bool ExtraOK;
	for (unsigned int l=0; l<(l.Extras).size(); l++) {
		ExtraOK = false;
		for (unsigned int m=0; m<(l.Extras).size(); m++) {
			if (((l.Extras[m]).first).compare(((u.Extras[l]).first))==0) ExtraOK = true;
		}
		if (!ExtraOK) return false;
	}
	return true;
}

// constructor(lower,upper)
FilenameComparator::FilenameComparator(const FilenameAttributes& l, const FilenameAttributes& u): Lower(l), Upper(u) {
	check(l,u);
}

// operator=
FilenameComparator& FilenameComparator::operator=(const FilenameComparator& rhs) {
	copy(rhs);
	return *this;
}

// set
void set(const FilenameAttributes& l, const FilenameAttributes& u) {
	if (check(l,u)) {
		Lower = l;
		Upper = u;
	}
}

// setLower
void setLower(const FilenameAttributes& l) {
	if (check(l,Upper)) {
		Lower = l;
	}
}

// setUpper
void setUpper(const FilenameAttributes& u) {
	if (check(Lower,u)) {
		Upper = u;
	}
}

// operator(Filename)
bool FilenameComparator::operator()(const Filename& f) const{
	if (!(Lower.Directory).empty()) {
		if ((f.Directory).compare(Lower.Directory)!=0) return false;
	}
	if (!(Lower.Timenumber).empty()) {
		if (f.Timenumber<Lower.Timenumber || f.Timenumber > Upper.Timenumber) return false;
	}
	if (!(Lower.ID).empty()) {
		if ((f.ID).compare(Lower.ID)!=0) return false;
	}
	if (!(Lower.Suffix).empty()) {
		if ((f.Suffix).compare(Lower.Suffix)!=0) return false;
	}
	size_t NumExtras = (Lower.Extras).size();
	if (NumExtras>0) {
		if ((f.Extras).size()!=NumExtras) return false;
		bool ExtraOK;
		for (unsigned int l=0; l<NumExtras; l++) {
			ExtraOK = false;
			for (unsigned int m=0; m<NumExtras; m++) {
				if (((Lower.Extras[m]).first).compare(((f.Extras[l]).first))==0) {
					if (((f.Extras[l]).second)>=((Lower.Extras[m]).second) && ((f.Extras[l]).second)<=((Upper.Extras[m]).second))
						ExtraOK = true;
				}
			}
			if (!ExtraOK) return false;
		}
	}
	return true;
}

// <<
ostream& operator<<(ostream& os, const FilenameComparator& fc){
	os << "Lower: " << endl << fc.Lower << endl << "Upper: " << fc.Upper << endl;
}

/*
	definitions for the Folder class.
		- update
		- <<
*/

void Folder::update() {
	int systemCall = system("dir data/* > dataFiles");
	if (systemCall==-1)
		cerr << "system call failure, finding dataFiles" << endl;
}

