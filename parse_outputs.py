#!/usr/bin/env python3
import argparse
import sys

def main():
    parser = argparse.ArgumentParser( description = 'Parse output produced by netMHCpan-4.0' )

    parser.add_argument( '-f', '--in_file',
                         help = "Output file produced by netMHCpan-4.0 to parse"
                       )

    parser.add_argument( '-w', '--weak',
                         help = (
                                  "Floating-point threshold for a binding to be considered "
                                  "weak. Anything less than this value will be considered a weak "
                                  "binding if it is also greater than or equal "
                                  "the strong threshold. (0.5, 0.2] "
                                )
                       )
    parser.add_argument( '-s', '--strong',
                         help = (
                                  "Floating-point threshold for a binding to be considered "
                                  "strong. Anything less than this, "
                                  "will be considered a strong binding. "
                                  "(0, 0.2)"
                                )
                       )
    parser.add_argument( '-o', '--output', default = "output.txt",
                         help = "File to write tab-delimited output to"
                       )
    parser.add_argument( '-v', '--verbose', default = False,
                         action = "store_true",
                         help = "Flag for output to be written to STDOUT as well as output file"
                       )
                         

    args = parser.parse_args()

    HEADER = "Pep-len\tAllele\tTotal\t#SB\t#WB\n" 

    output_parser = BindingOutputParser( args.in_file )
    output_parser.read_file()

    if not output_parser.file_found():
        print( "ERROR: %s was either unable to be found or opened, exiting" % output_parser.get_file_name() )
        sys.exit( 1 )
    allele_info = output_parser.parse()

    for current_allele in allele_info:
        current_allele.set_weak_binding_thresh( float( args.weak ) )
        current_allele.set_strong_binding_thresh( float( args.strong ) )
    

    out_file = open( args.output, 'w' )
    out_file.write( HEADER )

    if args.verbose:
        # Don't print a newline 
        print( HEADER[ :-1: ] )

    for current_allele in allele_info:
        out_file.write( str( current_allele ) + '\n'  )
        if args.verbose:
            print( str( current_allele ) )

    out_file.close()

class AlleleInfo:
    def __init__( self, name, weak_binding_thresh = 0, strong_binding_thresh = 0 ):
        self._name = name
        self._weak_binding_thresh = weak_binding_thresh
        self._strong_binding_thresh = strong_binding_thresh
        self._bindings = list()
        self._seq_names = set()
        self._lengths = set()
        self._header = "Pep-len\tAllele\tTotal\t#SB\t#WB"


    class BindingInfo:
        def __init__( self, attributes_list ):
            self._pos = attributes_list[ 0 ]
            self._peptide = attributes_list[ 1 ]
            self._peptide_id = attributes_list[ 2 ]
            self._core = attributes_list[ 3 ]
            self._icore = attributes_list[ 4 ]
            self._log = attributes_list[ 5 ]
            self._nm = attributes_list[ 6 ]
            self._rank = attributes_list[ 7 ]
            self._peptide_length = len( attributes_list[ 1 ] )
            self._attributes = attributes_list

        def is_strong( self, strong_threshold, weak_threshold ):
            return float( self._rank ) < strong_threshold

        def is_weak( self, strong_threshold, weak_threshold ):
            return float( self._rank ) <= weak_threshold and float( self._rank ) > strong_threshold

        def __str__( self ):
            out_str = "\t".join( self._attributes )
            return out_str
    def set_weak_binding_thresh( self, new_thresh ):
        self._weak_binding_thresh = new_thresh

    def set_strong_binding_thresh( self, new_thresh ):
        self._strong_binding_thresh = new_thresh

    def get_name( self ):
        return self._name

    def get_num_bindings_of_length( self, length ):
        total = 0
        for current_binding in self._bindings:
            if current_binding._peptide_length == length:
                total += 1
        return total

    def get_strong_bindings( self ):
        out_list = list()
        for current_binding in self._bindings:
            if current_binding.is_strong( self._strong_binding_thresh, self._weak_binding_thresh ):
                out_list.append( current_binding )
        return out_list

    def get_strong_bindings_by_length( self, length ):
        strong_bindings = self.get_strong_bindings()
        out_list = list()

        for current_binding in strong_bindings:
            if current_binding._peptide_length == length:
                out_list.append( current_binding )
        return out_list

    def get_weak_bindings_by_length( self, length ):
       weak_bindings = self.get_weak_bindings()
       out_list = list()

       for current_binding in weak_bindings:
           if current_binding._peptide_length == length:
               out_list.append( current_binding )
       return out_list               
        
    def get_weak_bindings( self ):
       out_list = list()
       for current_binding in self._bindings:
           if current_binding.is_weak( self._strong_binding_thresh, self._weak_binding_thresh ):
               out_list.append( current_binding )
       return out_list
           
    def get_strong_bindings_per_seq( self ):
        strong_bindings = self.get_strong_bindings()
        out_bindings = {}
        for current_binding in strong_bindings:
            if current_binding not in out_bindings:
                out_bindings[ current_binding._peptide_id ] = 0
            out_bindings[ current_binding._peptide_id ] += 1
        return out_bindings
                
    def get_weak_bindings_per_seq( self ):
       weak_bindings = self.get_weak_bindings()
       out_bindings = {}
       for current_binding in weak_bindings:
           if current_binding not in out_bindings:
               out_bindings[ current_binding._peptide_id ] = 0
           out_bindings[ current_binding._peptide_id ] += 1
       return out_bindings           
            


    @staticmethod
    def _is_header( in_list ):
        return in_list[ 0 ] == "Pos"

    def add_info( self, string ):
        split_string = string.split()
        if not AlleleInfo._is_header( split_string ):
            binding_info = self.BindingInfo( split_string )

            self._seq_names.add( binding_info._peptide_id )
            self._bindings.append( binding_info )
            self._lengths.add( binding_info._peptide_length )
        else:
            self._header = "\t".join( split_string )

    def __str__( self ):
       output_string = ""
       sorted_lengths = sorted( self._lengths )
       for current_length in sorted_lengths:
           num_weak_bindings = len( self.get_weak_bindings_by_length( current_length ) )
           num_strong_bindings = len( self.get_strong_bindings_by_length( current_length ) )
           total_bindings = num_weak_bindings + num_strong_bindings

           output_string += "%d\t%s\t%d\t%s\t%s\n" % (  current_length,
                                                    self._name,
                                                    total_bindings,
                                                    num_strong_bindings,
                                                    num_weak_bindings,
                                                 )
           
       # No newline at end of string 
       return output_string[:-1:]

class BindingOutputParser:
    def __init__( self, in_file ):
        self._file_name = in_file
        self._lines_from_file = None
        self._in_file_not_found = False
        self._alleles_from_file = list()

    def read_file( self ):
        file_opened = False
        try:
            open_file = open( self._file_name, 'r' )
            self._lines_from_file = open_file.readlines()
            sucess = True
        except ( IOError, OSError, TypeError ):
            self._in_file_not_found = True
        return file_opened

    def get_file_name( self ):
        return self._file_name
        
    def file_found( self ):
       return not self._in_file_not_found       
        
    def get_file_name( self ):
        return self._file_name

    def get_alleles( self ):
        output_alleles = ""
        if self._lines_from_file is not None:
            output_alleles = self._lines_from_file[ 0 ]
            output_alleles = output_alleles.strip().split()
            
        return output_alleles

    def add_allele( self, to_add):
        allele_to_add = AlleleInfo( to_add )
        self._alleles_from_file.append( allele_to_add )

    def create_alleles( self ):
        allele_names = self.get_alleles()

        if len( allele_names ) > 0:
            for current in allele_names:
                self.add_allele( current )
        
    
    def _parse_line( self, line_to_parse ):
        line_split = line_to_parse.split()
        output = list()
        num_alleles = len( self._alleles_from_file )

        for current in range( num_alleles ):
            current_list = list()
            current_list.append( line_split[ 0 ] )
            current_list.append( line_split[ 1 ] )
            current_list.append( line_split[ 2 ] )

            current_list.append( line_split[ ( current * 5 ) + 3 ] )
            current_list.append( line_split[ ( current * 5 ) + 4 ] )
            current_list.append( line_split[ ( current * 5 ) + 5 ] )
            current_list.append( line_split[ ( current * 5 ) + 6 ] )
            current_list.append( line_split[ ( current * 5 ) + 7 ] )

            output.append( current_list )
        return output



    def parse( self ):
        allele_info_list = list()
        if self.file_found():
            self.create_alleles()
            allele_info_list = self._alleles_from_file
            num_alleles = len( self._alleles_from_file )

            for current_line in self._lines_from_file[ 2:: ]:
                line_stripped = current_line.strip()
                parsed_line = self._parse_line( current_line )

                for current in range( len( parsed_line ) ):
                    allele_info_list[ current ].add_info( '\t'.join( parsed_line[ current ] ) )
                    
        return allele_info_list

    def get_info( self ):
        self.read_file()
        return self.parse()

if __name__ == '__main__':
    main()

