CC     := cc
CFLAG  := -std=c99 -Wall -Wextra -O3
INC    := -Iinclude
LIB    := -lm
SRCDIR := src
OBJDIR := obj
SRCS   := $(shell find $(SRCDIR) -type f -name *.c)
OBJS   := $(subst $(SRCDIR)/,$(OBJDIR)/,$(SRCS:.c=.o))
DEPS   := $(subst $(SRCDIR)/,$(OBJDIR)/,$(SRCS:.c=.d))
OUTDIR := output
TARGET := a.out

help:
	@echo "all     : create \"$(TARGET)\" with implicit scheme without X"
	@echo "all-x   : create \"$(TARGET)\" with implicit scheme with X"
	@echo "exp     : create \"$(TARGET)\" with explicit scheme without X"
	@echo "exp-x   : create \"$(TARGET)\" with explicit scheme with X"
	@echo "clean   : remove \"$(TARGET)\" and object files under \"$(OBJDIR)\""
	@echo "output  : create \"$(OUTDIR)\" to store output"
	@echo "datadel : clean-up \"$(OUTDIR)\""
	@echo "help    : show this message"

all: CFLAG += -DSCHEME=0
all: $(TARGET)

exp: CFLAG += -DSCHEME=1
exp: $(TARGET)

all-x: CFLAG += -DENABLE_XWINDOW
all-x: LIB   += -lX11
all-x: $(TARGET)

exp-x: CFLAG += -DSCHEME=1
exp-x: CFLAG += -DENABLE_XWINDOW
exp-x: LIB   += -lX11
exp-x: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAG) -o $@ $^ $(LIB)

$(OBJDIR)/%.o: $(SRCDIR)/%.c
	@if [ ! -e `dirname $@` ]; then \
		mkdir -p `dirname $@`; \
	fi
	$(CC) $(CFLAG) -MMD $(INC) -c $< -o $@

clean:
	$(RM) -r $(OBJDIR) $(TARGET)

output:
	@if [ ! -e $(OUTDIR) ]; then \
	   mkdir $(OUTDIR); \
	fi

datadel:
	@$(RM) $(OUTDIR)/*.dat

-include $(DEPS)

.PHONY : all clean output datadel help

