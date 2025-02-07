(define-module (celegans_genotype_converter)
	       #:use-module (guix packages) 
	       #:use-module (guix git-download) 
	       #:use-module (guix build-system python) 
	       #:use-module (guix licenses) 
	       #:use-module (gnu packages python) 
	       #:use-module (gnu packages base)) 

(define-public celegans_genotype_converter
	       (package 
		 (name "celegans_genotype_converter")
		 (version "0.1.0") 
		 (source 
		   (origin
		     (method git-fetch) 
		     (uri (git-reference 
			    (url "https://github.com/fetche-lab/GeneNetwork_24L/tree/main/Data_processing/C_elegans/Scripts/celegans_geno_converter.git")
			    (commit "3d613e6a4b4aab48e618f3e6525eb22ca03e19a7")))
		     (sha256
		       (base32 "sha256-1a29cs53qi3ya3icln5dwzmhm56crhcpy802s2kifg6sf88zqk2j"))))
		 (build-system python-build-system)
		 (inputs
		   `(("python" ,python-3)))
		 (synopsis "A package for genotype data conversion") 
		 (description " A Python package that converts genotype data into GEMMA-compatible format.")
		 (home-page "https://github.com/fetche-lab/GeneNetwork_24L/tree/main/Data_processing/C_elegans/Scripts/celegans_geno_converter.git")
		 (license mit)))
