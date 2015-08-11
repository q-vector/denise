//
// symbolic.h
// 
// Copyright (C) 2005-2013 Simon E. Ching
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

#ifndef DENISE_SYMBOLIC_H
#define DENISE_SYMBOLIC_H

#include <cmath>
#include <stack>
#include <denise/dstring.h>

using namespace std;
using namespace denise;

/// denise
namespace denise
{

   class Expression
   {

      public:

         class Instance
         {

            public:

               virtual Real
               evaluate (const Dstring& variable) const = 0;

         };

      private:

         class Node
         {

            public:

               virtual Real
               evaluate (const Instance& instance) const = 0;

         };

         class Constant : public Node
         {

            private:

               Real
               value;

            public:

               Constant (const Dstring& token);

               Real
               evaluate (const Instance& instance) const;

         };

         class Variable : public Node,
                          public Dstring
         {

            public:

               Variable (const Dstring& token);

               Real
               evaluate (const Instance& instance) const;

         };

         class Mono_Op : public Node
         {

            protected:

               Node*
               node_ptr;

            public:

               Mono_Op (Node* node_ptr);

               static Mono_Op*
               create (const Dstring& token,
                       Node* node_ptr);

         };

         class Bi_Op : public Node
         {

            protected:

               Node*
               left_ptr;

               Node*
               right_ptr;

            public:

               Bi_Op (Node* left_ptr,
                      Node* right_ptr);

               static Bi_Op*
               create (const Dstring& token,
                       Node* left_ptr,
                       Node* right_ptr);

         };

         class Exp : public Mono_Op
         {

            public:

               Exp (Node* node_ptr);

               Real
               evaluate (const Instance& instance) const;

         };

         class Log : public Mono_Op
         {

            public:

               Log (Node* node_ptr);

               Real
               evaluate (const Instance& instance) const;

         };

         class Sin : public Mono_Op
         {

            public:

               Sin (Node* node_ptr);

               Real
               evaluate (const Instance& instance) const;

         };

         class Cos : public Mono_Op
         {

            public:

               Cos (Node* node_ptr);

               Real
               evaluate (const Instance& instance) const;

         };

         class Tan : public Mono_Op
         {

            public:

               Tan (Node* node_ptr);

               Real
               evaluate (const Instance& instance) const;

         };

         class Sinh : public Mono_Op
         {

            public:

               Sinh (Node* node_ptr);

               Real
               evaluate (const Instance& instance) const;

         };

         class Cosh : public Mono_Op
         {

            public:

               Cosh (Node* node_ptr);

               Real
               evaluate (const Instance& instance) const;

         };

         class Tanh : public Mono_Op
         {

            public:

               Tanh (Node* node_ptr);

               Real
               evaluate (const Instance& instance) const;

         };

         class Asin : public Mono_Op
         {

            public:

               Asin (Node* node_ptr);

               Real
               evaluate (const Instance& instance) const;

         };

         class Acos : public Mono_Op
         {

            public:

               Acos (Node* node_ptr);

               Real
               evaluate (const Instance& instance) const;

         };

         class Atan : public Mono_Op
         {

            public:

               Atan (Node* node_ptr);

               Real
               evaluate (const Instance& instance) const;

         };

         class Asinh : public Mono_Op
         {

            public:

               Asinh (Node* node_ptr);

               Real
               evaluate (const Instance& instance) const;

         };

         class Acosh : public Mono_Op
         {

            public:

               Acosh (Node* node_ptr);

               Real
               evaluate (const Instance& instance) const;

         };

         class Atanh : public Mono_Op
         {

            public:

               Atanh (Node* node_ptr);

               Real
               evaluate (const Instance& instance) const;

         };

         class Plus : public Bi_Op
         {

            public:

               Plus (Node* left_ptr,
                     Node* right_ptr);

               Real
               evaluate (const Instance& instance) const;

         };

         class Minus : public Bi_Op
         {

            public:

               Minus (Node* left_ptr,
                      Node* right_ptr);

               Real
               evaluate (const Instance& instance) const;

         };

         class Multiply : public Bi_Op
         {

            public:

               Multiply (Node* left_ptr,
                         Node* right_ptr);

               Real
               evaluate (const Instance& instance) const;

         };

         class Divide : public Bi_Op
         {

            public:

               Divide (Node* left_ptr,
                       Node* right_ptr);

               Real
               evaluate (const Instance& instance) const;

         };

         class Power : public Bi_Op
         {

            public:

               Power (Node* left_ptr,
                      Node* right_ptr);

               Real
               evaluate (const Instance& instance) const;

         };

         std::stack<Node*>
         node_ptr_stack;

         set<Dstring>
         variable_set;

         static bool
         is_constant (const Dstring& token);

         static bool
         is_op (const Dstring& token);

      public:

         Expression (const Tokens& postfix_tokens);

   };

}

#endif /* DENISE_SYMBOLIC_H */
