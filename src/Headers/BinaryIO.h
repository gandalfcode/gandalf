//=================================================================================================
//  BinaryIO.h
//  Helper classes to write and read from binary files
//
//  This file is part of GANDALF :
//  Graphical Astrophysics code for N-body Dynamics And Lagrangian Fluids
//  https://github.com/gandalfcode/gandalf
//  Contact : gandalfcode@gmail.com
//
//  Copyright (C) 2013  D. A. Hubber, G. Rosotti
//
//  GANDALF is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 2 of the License, or
//  (at your option) any later version.
//
//  GANDALF is distributed in the hope that it will be useful, but
//  WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  General Public License (http://www.gnu.org/licenses) for more details.
//=================================================================================================


//=================================================================================================
//  Class BinaryWriter
/// \brief  Helper class to write to binary files
/// \author G. Rosotti
/// \date   13/02/2014
//=================================================================================================
class BinaryWriter {
private:
  std::ofstream& output_file;

public:
  BinaryWriter(std::ofstream& _output_file) :
    output_file(_output_file) {}

  void write_value(const string& value) {
    output_file << value;
  }

  void write_value(const bool& value) {
    const int integer_value = value;
    write_value(integer_value);
  }

  template <class TypeToWrite>
  void write_value(const TypeToWrite& value) {
    output_file.write((char*) &value, sizeof(value));
  }

};



//=================================================================================================
//  Class BinaryReader
/// \brief  Helper class for reading from binary files
/// \author G. Rosotti
/// \date   13/02/2014
//=================================================================================================
class BinaryReader {
private:
  std::ifstream& input_file;

public:
  BinaryReader(std::ifstream& _input_file) :
    input_file (_input_file) {}

  template <class TypeToRead>
  void read_value(TypeToRead& result) {
    input_file.read((char*) &result, sizeof(TypeToRead));
  }

};
