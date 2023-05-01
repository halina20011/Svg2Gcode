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

typedef struct Node{
    void *data;
    struct Node *next;
} Node;

int initSingleLinkedList(Node **head, Node **tail){
    *head = malloc(sizeof(Node));
    if(*head == NULL){
        return 1;
    }

    (*head)->data = NULL;
    (*head)->next = NULL;
    *tail = *head;

    return 0;
}

void printSingleLinkedList(Node *head){
    Node *current = head;
    while((current = current->next) != NULL){
        char *s = (void*)(current->data);
        printf("%s\n", s);
    }
}

// Add node to the end of the list
Node *add(Node *head, void *data){
    Node *newNode = malloc(sizeof(Node));
    if(newNode == NULL){
        return NULL;
    }

    newNode->data = data;
    newNode->next = NULL;

    Node *current = head;
    while(current->next != NULL){
        current = current->next;
    }

    current->next = newNode;

    return newNode;
}

// if both are null: then return null
// if tail is null: get last node from head
// if head->next null: set it to made node
Node *push(Node *head, Node **tail, void *data){
    if(head == NULL || tail == NULL){
        return NULL;
    }

    Node *newNode = malloc(sizeof(Node));
    if(newNode == NULL){
        return NULL;
    }

    newNode->data = data;
    newNode->next = NULL;

    if(*tail == NULL){
        *tail = head;
        while((*tail = (*tail)->next) != NULL){
        }
    }

    (*tail)->next = newNode;
    *tail = newNode;

    if(head->next == NULL){
        printf("Head n is NULL\n");
        head->next = newNode;
    }

    return newNode;
}

Node *insert(Node *head, char *data){
    Node *newNode = malloc(sizeof(Node));
    if(newNode == NULL){
        return NULL;
    }

    newNode->next = head->next;
    newNode->data = data;

    head->next = newNode;

    return newNode;
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

void freeSingleLinkedList(Node *head){
    Node *current = head->next;

    while(current != NULL){
        Node *next = current->next;
        free(current->data);
        free(current);
        current = next;
    }
    
    head->next = NULL;

    free(head);
}

#endif
