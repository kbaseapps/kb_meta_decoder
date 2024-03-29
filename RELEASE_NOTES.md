# kb_meta_decoder release notes
=========================================

1.0.2
-----
* Fixed bug in which MIN_DEPTH was treated as MAX_DEPTH
* Added new MAX_DEPTH parameter to give users more control over DP filtering

1.0.1
-----
* Updated documentation
* Made "calculate population statistics" inactive until upstream is fixed

1.0.0
-----
* Updated figure generation for "call snps"
* Added more documentation
* Initial non-beta release

0.0.6
-----
* Updated reports; limited debug info in Summary
* Added new icons
* Changed name of "Call variants" method to "Call SNPs"

0.0.5
-----
* Moved each job to separate output dir, to avoid filename collisions in scratch

0.0.4
-----
* Fixed unit tests to work with blob store, and added several small tests

0.0.3
-----
* Call variants does multiple reads libraries at a time, in parallel

0.0.2
-----
* Initial version for public beta testing

0.0.1
-----
* Initial dev version

0.0.0
-----
* Module created by kb-sdk init
