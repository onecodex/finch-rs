# Copyright © 2015, Battelle National Biodefense Institute (BNBI);
# all rights reserved. Authored by: Brian Ondov, Todd Treangen,
# Sergey Koren, and Adam Phillippy
#
# See the LICENSE.txt file included with this software for license information.
# Content of LICENSE.txt
#PURPOSE
#
#Mash is a fast sequence distance estimator that uses the MinHash
#algorithm and is designed to work with genomes and metagenomes in the
#form of assemblies or reads. It is implemented in C++ and is
#distributed with:
#
#KSeq
#  lh3lh3.users.sourceforge.net/kseq.shtml
#  MIT License
#
#MurmurHash3
#  code.google.com/p/smhasher/wiki/MurmurHash3
#  Public domain
#
#Open Bloom Filter
#  https://code.google.com/p/bloom/source/browse/trunk/bloom_filter.hpp
#  Common Public License
#
#COPYRIGHT LICENSE
#
#Copyright © 2015, Battelle National Biodefense Institute (BNBI);
#all rights reserved. Authored by: Brian Ondov, Todd Treangen,
#Sergey Koren, and Adam Phillippy
#
#This Software was prepared for the Department of Homeland Security
#(DHS) by the Battelle National Biodefense Institute, LLC (BNBI) as
#part of contract HSHQDC-07-C-00020 to manage and operate the National
#Biodefense Analysis and Countermeasures Center (NBACC), a Federally
#Funded Research and Development Center.
#
#Redistribution and use in source and binary forms, with or without
#modification, are permitted provided that the following conditions are
#met:
#
#1. Redistributions of source code must retain the above copyright
#notice, this list of conditions and the following disclaimer.
#
#2. Redistributions in binary form must reproduce the above copyright
#notice, this list of conditions and the following disclaimer in the
#documentation and/or other materials provided with the distribution.
#
#3. Neither the name of the copyright holder nor the names of its
#contributors may be used to endorse or promote products derived from
#this software without specific prior written permission.
#
#THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
#"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
#LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
#A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
#HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
#SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
#LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
#DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
#THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
#(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
#OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# + https://github.com/marbl/Mash/issues/112

using Cxx = import "/capnp/c++.capnp";
$Cxx.namespace("capnp");

@0xc4c8b1ada05e7704;

struct MinHash
{
	struct ReferenceList
	{
		struct Reference
		{
			sequence @0 : Text;
			quality @1 : Text;
			length @2 : UInt32;
			length64 @ 7 : UInt64;
			# https://github.com/marbl/Mash/issues/112
			numValidKmers @9 : UInt64;
			name @3 : Text;
			comment @4 : Text;
			hashes32 @5 : List(UInt32);
			hashes64 @6 : List(UInt64);
			counts32 @8 : List(UInt32);
		}

		references @0 : List(Reference);
	}

	struct LocusList
	{
		struct Locus
		{
			sequence @0 : UInt32;
			position @1 : UInt32;
			hash32 @2 : UInt32;
			hash64 @3 : UInt64;
		}

		loci @0 : List(Locus);
	}

	kmerSize @0 : UInt32;
	windowSize @1 : UInt32;
	minHashesPerWindow @2 : UInt32;
	concatenated @3 : Bool;
	error @6 : Float32;
	noncanonical @7 : Bool;
	alphabet @8 : Text;
	preserveCase @9 : Bool;
	hashSeed @10 : UInt32 = 42;

	referenceListOld @4 : ReferenceList;
	referenceList @11 : ReferenceList;
	locusList @5 : LocusList;
}
