// symbolic.cc
// 
// Copyright (C) 2010 Simon E. Ching
// 
// This file is part of libdenise.
//
// libdenise is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// libdenise is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with libdenise.  If not, see <http://www.gnu.org/licenses/>.

#include "symbolic.h"

Expression::Constant::Constant (const wstring& token)
   : value (stof (token))
{
}

Real
Expression::Constant::evaluate (const Instance& instance) const
{
   return value;
}

Expression::Variable::Variable (const wstring& token)
   : wstring (token)
{
}

Real
Expression::Variable::evaluate (const Instance& instance) const
{
   return instance.evaluate (*this);
}

Expression::Mono_Op::Mono_Op (Node* node_ptr)
   : node_ptr (node_ptr)
{
}

Expression::Mono_Op*
Expression::Mono_Op::create (const wstring& token,
                             Node* node_ptr)
{
   if (token == L"exp") { return new Exp (node_ptr); }
   if (token == L"log") { return new Log (node_ptr); }
   if (token == L"sin") { return new Sin (node_ptr); }
   if (token == L"cos") { return new Cos (node_ptr); }
   if (token == L"tan") { return new Tan (node_ptr); }
   if (token == L"sinh") { return new Sinh (node_ptr); }
   if (token == L"cosh") { return new Cosh (node_ptr); }
   if (token == L"tanh") { return new Tanh (node_ptr); }
   if (token == L"asin") { return new Asin (node_ptr); }
   if (token == L"acos") { return new Acos (node_ptr); }
   if (token == L"atan") { return new Atan (node_ptr); }
   if (token == L"asinh") { return new Asinh (node_ptr); }
   if (token == L"acosh") { return new Acosh (node_ptr); }
   if (token == L"atanh") { return new Atanh (node_ptr); }
   throw Exception (L"Not Mono_Op");
}

Expression::Bi_Op::Bi_Op (Node* left_ptr,
                          Node* right_ptr)
   : left_ptr (left_ptr),
     right_ptr (right_ptr)
{
}

Expression::Bi_Op*
Expression::Bi_Op::create (const wstring& token,
                           Node* left_ptr,
                           Node* right_ptr)
{
   if (token == L"+") { return new Plus (left_ptr, right_ptr); }
   if (token == L"-") { return new Minus (left_ptr, right_ptr); }
   if (token == L"*") { return new Multiply (left_ptr, right_ptr); }
   if (token == L"/") { return new Divide (left_ptr, right_ptr); }
   if (token == L"^") { return new Power (left_ptr, right_ptr); }
   throw Exception (L"Not Bi_Op");
}

Expression::Exp::Exp (Node* node_ptr)
   : Mono_Op (node_ptr)
{
}
 
Real
Expression::Exp::evaluate (const Instance& instance) const
{
   return exp (node_ptr->evaluate (instance));
}

Expression::Log::Log (Node* node_ptr)
   : Mono_Op (node_ptr)
{
}
 
Real
Expression::Log::evaluate (const Instance& instance) const
{
   return log (node_ptr->evaluate (instance));
}

Expression::Sin::Sin (Node* node_ptr)
   : Mono_Op (node_ptr)
{
}
 
Real
Expression::Sin::evaluate (const Instance& instance) const
{
   return sin (node_ptr->evaluate (instance));
}

Expression::Cos::Cos (Node* node_ptr)
   : Mono_Op (node_ptr)
{
}
 
Real
Expression::Cos::evaluate (const Instance& instance) const
{
   return cos (node_ptr->evaluate (instance));
}

Expression::Tan::Tan (Node* node_ptr)
   : Mono_Op (node_ptr)
{
}
 
Real
Expression::Tan::evaluate (const Instance& instance) const
{
   return tan (node_ptr->evaluate (instance));
}

Expression::Sinh::Sinh (Node* node_ptr)
   : Mono_Op (node_ptr)
{
}
 
Real
Expression::Sinh::evaluate (const Instance& instance) const
{
   return sinh (node_ptr->evaluate (instance));
}

Expression::Cosh::Cosh (Node* node_ptr)
   : Mono_Op (node_ptr)
{
}
 
Real
Expression::Cosh::evaluate (const Instance& instance) const
{
   return cosh (node_ptr->evaluate (instance));
}

Expression::Tanh::Tanh (Node* node_ptr)
   : Mono_Op (node_ptr)
{
}
 
Real
Expression::Tanh::evaluate (const Instance& instance) const
{
   return tanh (node_ptr->evaluate (instance));
}

Expression::Asin::Asin (Node* node_ptr)
   : Mono_Op (node_ptr)
{
}
 
Real
Expression::Asin::evaluate (const Instance& instance) const
{
   return asin (node_ptr->evaluate (instance));
}

Expression::Acos::Acos (Node* node_ptr)
   : Mono_Op (node_ptr)
{
}
 
Real
Expression::Acos::evaluate (const Instance& instance) const
{
   return acos (node_ptr->evaluate (instance));
}

Expression::Atan::Atan (Node* node_ptr)
   : Mono_Op (node_ptr)
{
}
 
Real
Expression::Atan::evaluate (const Instance& instance) const
{
   return atan (node_ptr->evaluate (instance));
}

Expression::Asinh::Asinh (Node* node_ptr)
   : Mono_Op (node_ptr)
{
}
 
Real
Expression::Asinh::evaluate (const Instance& instance) const
{
   return asinh (node_ptr->evaluate (instance));
}

Expression::Acosh::Acosh (Node* node_ptr)
   : Mono_Op (node_ptr)
{
}
 
Real
Expression::Acosh::evaluate (const Instance& instance) const
{
   return acosh (node_ptr->evaluate (instance));
}

Expression::Atanh::Atanh (Node* node_ptr)
   : Mono_Op (node_ptr)
{
}
 
Real
Expression::Atanh::evaluate (const Instance& instance) const
{
   return atanh (node_ptr->evaluate (instance));
}

Expression::Plus::Plus (Node* left_ptr,
                        Node* right_ptr)
   : Bi_Op (left_ptr, right_ptr)
{
}

Real
Expression::Plus::evaluate (const Instance& instance) const
{
   return left_ptr->evaluate (instance) + right_ptr->evaluate (instance);
}

Expression::Minus::Minus (Node* left_ptr,
                          Node* right_ptr)
   : Bi_Op (left_ptr, right_ptr)
{
}

Real
Expression::Minus::evaluate (const Instance& instance) const
{
   return left_ptr->evaluate (instance) - right_ptr->evaluate (instance);
}

Expression::Multiply::Multiply (Node* left_ptr,
                                Node* right_ptr)
   : Bi_Op (left_ptr, right_ptr)
{
}

Real
Expression::Multiply::evaluate (const Instance& instance) const
{
   return left_ptr->evaluate (instance) * right_ptr->evaluate (instance);
}

Expression::Divide::Divide (Node* left_ptr,
                            Node* right_ptr)
   : Bi_Op (left_ptr, right_ptr)
{
}

Real
Expression::Divide::evaluate (const Instance& instance) const
{
   return left_ptr->evaluate (instance) / right_ptr->evaluate (instance);
}

Expression::Power::Power (Node* left_ptr,
                          Node* right_ptr)
   : Bi_Op (left_ptr, right_ptr)
{
}

Real
Expression::Power::evaluate (const Instance& instance) const
{
   return pow (left_ptr->evaluate (instance), right_ptr->evaluate (instance));
}

bool
Expression::is_constant (const wstring& token)
{
   return Reg_Exp (L"[0-9]+").match (token);
}

bool
Expression::is_op (const wstring& token)
{
   if (token == L"+") { return true; }
   if (token == L"-") { return true; }
   if (token == L"*") { return true; }
   if (token == L"/") { return true; }
   if (token == L"^") { return true; }
   if (token == L"exp") { return true; }
   if (token == L"log") { return true; }
   if (token == L"sin") { return true; }
   if (token == L"cos") { return true; }
   if (token == L"tan") { return true; }
   if (token == L"sinh") { return true; }
   if (token == L"cosh") { return true; }
   if (token == L"tanh") { return true; }
   if (token == L"asin") { return true; }
   if (token == L"acos") { return true; }
   if (token == L"atan") { return true; }
   if (token == L"asinh") { return true; }
   if (token == L"acosh") { return true; }
   if (token == L"atanh") { return true; }
   return false;
}

Expression::Expression (const Tokens& postfix_tokens)
{

   for (Tokens::const_iterator iterator = postfix_tokens.begin ();
        iterator != postfix_tokens.end (); iterator++)
   {

      const wstring& token = *(iterator);

      if (is_op (token))
      {

         const Integer n = node_ptr_stack.size ();

         if (node_ptr_stack.size () < 1) { throw Exception (L"Empty Stack"); }
         Node* node_ptr = node_ptr_stack.top ();
         node_ptr_stack.pop ();

         Node* new_node_ptr = Mono_Op::create (token, node_ptr);

         if (new_node_ptr != NULL)
         {
            node_ptr_stack.push (new_node_ptr);
            continue;
         }
         else
         {
            if (node_ptr_stack.size () < 1) { throw Exception (L"Empty Stack"); }
            Node* left_ptr = node_ptr_stack.top ();
            node_ptr_stack.pop ();
            Node* new_node_ptr = Bi_Op::create (token, left_ptr, node_ptr);
            if (new_node_ptr == NULL) { throw Exception (L"Unrecognized Op"); }
            node_ptr_stack.push (new_node_ptr);
         }

      }
      else
      {
         if (is_constant (token))
         {
            // is a constant
            Node* node_ptr = new Constant (token);
            node_ptr_stack.push (node_ptr);
         }
         else
         {
            // is a variable
            Node* node_ptr = new Variable (token);
            node_ptr_stack.push (node_ptr);
         }
      }

   }

}

