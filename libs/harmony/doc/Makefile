DOXYGEN=doxygen

DOCS = Users_Guide \
       Plugin_Guide \
       Internal_Guide
IMGS = images/feedback \
       images/session \
       images/multi-client \
       images/accept \
       images/reject \
       images/return \
       images/retry

.PHONY: all clean distclean

all: $(addsuffix .pdf, $(DOCS))

$(addsuffix .pdf, $(DOCS)): %.pdf: doxygen/Doxyfile.% doxygen/Doxyfile.base doxygen/%.dox doxygen/harmony.sty $(addsuffix .pdf, $(IMGS)) images/hserver1.png images/hserver2.png
	@echo Making the $*
	@$(DOXYGEN) $< && \
	./doxygen/postprocess.sh $* && \
	$(MAKE) -C $*.latex > /dev/null 2>&1 && \
	cp $*.latex/refman.pdf $@

$(addsuffix .pdf, $(IMGS)): %.pdf: %.svg
	inkscape --export-pdf=$@ $<

clean:
	rm -rf $(addsuffix .pdf, $(DOCS)) \
		$(addsuffix .html, $(DOCS)) \
		$(addsuffix .latex, $(DOCS)) \
		$(addsuffix .pdf, $(IMGS))

distclean: clean
	rm -rf *~ */*~ *.d

#
# Auto dependency creation
#
$(addsuffix .d, $(DOCS)): %.d: doxygen/Doxyfile.%
	@rm -f $@ && \
	echo -n $*.pdf: >> $@ && \
	sed ':a;N;$$!ba;s/\\\n/ /g' $< | grep INPUT | sed 's/INPUT.*=//' >> $@

-include $(addsuffix .d, $(DOCS))
