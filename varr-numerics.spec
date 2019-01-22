# vim:set ft=spec:

Name:     varr-numerics
Version:  0.1.0
Release:  1%{?dist}
Summary:  Experimental variable resolution primitive numerics
Source0:  %{name}-%{version}.tar.xz
License:  GPLv3
Group: Development/Libraries

BuildRequires: cmake3
BuildRequires: cmake3-data
BuildRequires: pkgconfig
BuildRequires: make
BuildRequires: gcc
BuildRequires: gsl-devel
Requires:      gsl

%description
%{summary}

%package devel
Summary: %{summary} -- development package
Requires: %{name} = %{version}
Requires: gsl-devel
%description devel
%{summary}

%prep
%setup -c -T -D -a 0 -n %{name}-%{version}

%build
%cmake3 %{name}-%{version}

%install
%make_install

%files
%doc %{name}-%{version}/README.md
%license %{name}-%{version}/LICENSE
%{_libdir}/*.so.*

%files devel
%doc %{name}-%{version}/README.md
%license %{name}-%{version}/LICENSE
%{_libdir}/*.so
%{_includedir}/*.h
%{_libdir}/pkgconfig/*.pc

%changelog
* Tue Jan 08 2019 Duncan Macleod <duncan.macleod@ligo.org> - 0.1.0-1
- First cut at packaging
