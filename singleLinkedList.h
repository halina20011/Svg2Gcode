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

#include <stdlib.h>

struct Node{
    void *data;
    struct Node *next;
};

int initSingleLinkedList(struct Node **head, struct Node **tail){
    *head = malloc(sizeof(struct Node));
    if(*head == NULL){
        return 1;
    }

    (*head)->data = NULL;
    (*head)->next = NULL;

    (*tail)->data = NULL;
    (*tail)->next = NULL;

    return 0;
}

void printSingleLinkedList(struct Node *head){
    struct Node *current = head;
    while((current = current->next) != NULL){
        char *s = (void*)(current->data);
        printf("%s\n", s);
    }
}

// if both are null: then return null
// if tail is null: get last node from head
// if head->next null: set it to made node
void push(struct Node **head, struct Node **tail, void *data){
    struct Node *newNode = malloc(sizeof(struct Node));
    if(newNode == NULL){
        return;
    }

    newNode->data = data;
    newNode->next = NULL;
    
    if(*head == NULL){
        *head = newNode;
        *tail = newNode;
    }
    else{
        (*tail)->next = newNode;
        *tail = newNode;
    }
}

void insert(struct Node **head, struct Node **tail, char *data){
    struct Node *newNode = malloc(sizeof(struct Node));
    if(newNode == NULL){
        return;
    }

    newNode->next = *head;
    newNode->data = data;

    if(*head == NULL){
        *head = newNode;
        *tail = newNode;
    }
    else{
        *head = newNode;
    }
}

char *delete(struct Node **head, int index){
    struct Node **indirect = head;
    char *data = NULL;

    int i = 0;
    for(; i < index && *indirect != NULL; i++){
        indirect = &((*indirect)->next);
    }
    
    if(i == index && *indirect != NULL){
        struct Node *temp = *indirect;
        data = temp->data;
        *indirect = temp->next;
        free(temp);
    }

    return data;
}

char *deleteFront(struct Node **head, struct Node **tail){
    if(*head == NULL){
        return NULL;
    }

    char *data = (*head)->data;
    struct Node *temp = *head;
    if(*head == *tail){
        *head = NULL;
        *tail = NULL;
    }
    else{
        *head = (*head)->next;
    }

    free(temp);

    return data;
}

void freeSingleLinkedList(struct Node *head){
    struct Node *curr = head;

    while(curr){
        struct Node *temp = curr;
        curr = curr->next;
        free(temp->data);
        free(temp);
    }
}

#endif
