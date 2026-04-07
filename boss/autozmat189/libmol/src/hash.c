#include <stdio.h>
#include <string.h>
#include <malloc.h>

#include "hash.h"

#define	HASHSIZE	4093
#define	STEPSIZE	1024

typedef struct _Bucket {
	int		n;
	struct _Bucket	*next;
	} Bucket;

static Bucket	*bucket[HASHSIZE];
static char	*strtable = NULL;
static int	cur_len = 0;
static int	max_len = 0;

static int	installstr(char *str);
static int	hash(char *str);

/* add a new string to a string table, "strtable" */
static int	installstr (char *str)
{
	register int	len, index;

	len = strlen(str) + 1;
	if (max_len - cur_len < len) {
		max_len += STEPSIZE;
		strtable = !strtable ? (char *)malloc(max_len) : (char *)realloc(strtable, max_len);
		if (!strtable) return -1;
	}
	strcpy(strtable + cur_len, str);
	index = cur_len;
	cur_len += len;
	return index;
}

/* The hash function hashes string into a bucket number */
static int	hash (char *str)
{
	register unsigned int	h;

	h = 0;
	while (*str) h = (h << 3) + *str++;
	return h % HASHSIZE;
}

void	CleanUpStrTable (void)
{
	int	i;
	Bucket	*next;

	free((char *)strtable);
	cur_len = max_len = 0;
	for(i=0;i<HASHSIZE;i++) {
		while (bucket[i]) {
			next = bucket[i]->next;
			free((char *)bucket[i]);
			bucket[i] = next;
		}
	}
}

/* search a string 'str' in the hash table */
int	LookUpStrTable (char *str)
{
	register Bucket	*p;
	register int	h;

	/* get hash value */
	h = hash(str);

	/* see if the string 'str' already exists in the hash table */
	for(p=bucket[h];p;p=p->next) {
		if (strcmp(strtable+p->n, str) == 0) return p->n;
	}

	/* install a new entry to the hash table */
	if (!(p = (Bucket *)malloc(sizeof(Bucket)))) return -1;
	p->n = installstr(str);
	p->next = bucket[h];
	bucket[h] = p;
	return p->n;
}

char	*GetStrValue (int n)
   	  	/* key */
{
	if (n < 0 || n >= cur_len) return (char *)NULL;
	return (strtable ? &strtable[n] : (char *)NULL);
}

