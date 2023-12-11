// Copyright (C) 2023  halina20011
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.

#ifndef SINGLELINKEDLIST
#define SINGLELINKEDLIST

struct Node{
    void *data;
    struct Node *next;
};

int initSingleLinkedList(struct Node **head, struct Node **tail);
void push(struct Node **head, struct Node **tail, void *data);
void printSingleLinkedList(struct Node *head);
void insert(struct Node **head, struct Node **tail, char *data);
void *deleteIndex(struct Node **head, int index);
char *deleteFront(struct Node **head, struct Node **tail);
void freeSingleLinkedList(struct Node *head);

#endif
