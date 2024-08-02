root_flags= `root-config --cflags --ldflags --libs --glibs`
compiler_flags=-std=c++11 -Wall -bind_at_load -Wl,-stack_size -Wl,4000000
includes= -I. -I../cmssw_code/ -I/usr/local/include -I /opt/homebrew/Cellar/boost/1.85.0/include
DEPS= template_utils.h ranlux.c ../cmssw_code/*
EXE= gen_zp_template gen_zp_template2d gen_xy_template gen_xy_template2d gen_xy_template3d test_template
ODIR=../bin
BOOST_FLAG=-DSI_PIXEL_TEMPLATE_USE_BOOST
LIBRARIES = -L/usr/local/lib

all: $(EXE)


gen_zp_template: gen_zp_template.cc $(DEPS)
	c++ $(BOOST_FLAG) $(compiler_flags) gen_zp_template.cc -o $(ODIR)/gen_zp_template  $(root_flags) $(includes) $(LIBRARIES)

gen_zp_template2d: gen_zp_template2d.cc $(DEPS)
	c++  $(BOOST_FLAG) $(compiler_flags) gen_zp_template2d.cc -o $(ODIR)/gen_zp_template2d  $(root_flags) $(includes)  $(LIBRARIES)

gen_xy_template: gen_xy_template.cc $(DEPS)
	c++ $(BOOST_FLAG)  $(compiler_flags) gen_xy_template.cc   -o $(ODIR)/gen_xy_template $(root_flags) $(includes)  $(LIBRARIES)

gen_xy_template2d: gen_xy_template2d.cc $(DEPS)
	c++ $(BOOST_FLAG)  $(compiler_flags) gen_xy_template2d.cc  -o $(ODIR)/gen_xy_template2d $(root_flags) $(includes)  $(LIBRARIES)
	
gen_xy_template3d: gen_xy_template3d.cc $(DEPS)
	c++ $(BOOST_FLAG)  $(compiler_flags) gen_xy_template3d.cc  -o $(ODIR)/gen_xy_template3d $(root_flags) $(includes)  $(LIBRARIES)

test_template: test_template.cc $(DEPS)
	c++ $(BOOST_FLAG)  $(compiler_flags) test_template.cc -o $(ODIR)/test_template $(root_flags) $(includes)  $(LIBRARIES)

load_template_test: load_template_test.cc $(DEPS)
	c++ load_template_test.cc -o $(ODIR)/load_template_test $(compiler_flags) $(root_flags) $(includes) $(BOOST_FLAG)

clean:
	rm -f $(ODIR)/*

.PHONY: all clean

