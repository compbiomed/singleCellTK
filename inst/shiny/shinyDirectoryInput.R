#########################
# SHINY DIRECTORY INPUT #
#########################
# source code: https://github.com/wleepang/shiny-directory-input
# LICENSE:
# BSD 3-Clause License
#
# Copyright (c) 2018, W. Lee Pang, PhD
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
#   * Redistributions of source code must retain the above copyright notice, this
# list of conditions and the following disclaimer.
#
# * Redistributions in binary form must reproduce the above copyright notice,
# this list of conditions and the following disclaimer in the documentation
# and/or other materials provided with the distribution.
#
# * Neither the name of the copyright holder nor the names of its
# contributors may be used to endorse or promote products derived from
# this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
#          SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


#' @name choose.dir
#'
#' @title Choose a Folder Interactively
#'
#' Display an OS-native folder selection dialog under Mac OS X, Linux GTK+ or
#' Windows.
#'
#' @param default which folder to show initially
#' @param caption the caption on the selection dialog
#' @param useNew boolean, selects the type of dialog shown in windows
#'
#' @details
#' Uses an Apple Script, Zenity or Windows Batch script to display an OS-native
#' folder selection dialog.
#'
#' For Apple Script, with \code{default = NA}, the initial folder selection
#' is determined by default behavior of the "choose folder" script. Otherwise,
#' paths are expanded with \code{\link{path.expand}}.
#'
#' For Linux, with \code{default = NA}, the initial folder selection is
#' determined by defaul behavior of the zenity script.
#'
#' The new windows batch script allows both initial folder and caption to be set.
#' In the old batch script for Windows the initial folder is always ignored.
#'
#' @return
#' A length one character vector, character NA if 'Cancel' was selected.
#'
#'
choose.dir = function(default = NA, caption = NA, useNew=TRUE) {
  if (Sys.info()['sysname'] == 'Darwin') {
    return(choose.dir.darwin(default = default, caption = caption))
  } else if (Sys.info()['sysname'] == 'Linux') {
    return(choose.dir.linux(default = default, caption = caption))
  } else if (Sys.info()['sysname'] == 'Windows') {
    # Use batch script to circumvent issue w/ `choose.dir`/`tcltk::tk_choose.dir`
    # window popping out unnoticed in the back of the current window
    return(choose.dir.windows(default = default, caption = caption, useNew = useNew))
  }
  return(paste("Error: don't know how to show a folder dialog in", Sys.info()['sysname']) )
}

#' @name choose.dir.darwin
#'
#' @title The apple version of the choose folder
#'
#' @seealso \code{\link{choose.dir}}
#'
#' @return
#' A length one character vector, character NA if 'Cancel' was selected.
#'
choose.dir.darwin <- function(default = NA, caption = NA) {
  command = 'osascript'
  args = '-e "POSIX path of (choose folder{{prompt}}{{default}})"'

  if (!is.null(caption) && !is.na(caption) && nzchar(caption)) {
    prompt = sprintf(' with prompt \\"%s\\"', caption)
  } else {
    prompt = ''
  }
  args = sub('{{prompt}}', prompt, args, fixed = TRUE)
  
  if (!is.null(default) && !is.na(default) && nzchar(default)) {
    default = sprintf(' default location \\"%s\\"', path.expand(default))
  } else {
    default = ''
  }
  args = sub('{{default}}', default, args, fixed = TRUE)
  
  suppressWarnings({
    path = system2(command, args = args, stderr = TRUE)
  })
  if (!is.null(attr(path, 'status')) && attr(path, 'status')) {
    # user canceled
    path = NA
  } else {
    # cut any extra output lines, like "Class FIFinderSyncExtensionHost ..."
    path = tail(path, n=1)
  }

  return(path)
}


#' @name choose.dir.linux
#'
#' @title The linux version of the choose folder
#'
#' @seealso \code{\link{choose.dir}}
#'
#' @return
#' A length one character vector, character NA if 'Cancel' was selected.
#'
choose.dir.linux <- function(default = NA, caption = NA) {
  command = 'zenity'
  args = '--file-selection --directory'

  if (!is.null(default) && !is.na(default) && nzchar(default)) {
    args = paste(args, sprintf('--filename="%s"', default))
  }

  if (!is.null(caption) && !is.na(caption) && nzchar(caption)) {
    args = paste(args, sprintf('--title="%s"', caption))
  }

  suppressWarnings({
    path = system2(command, args = args, stderr = TRUE)
  })

  #Return NA if user hits cancel
  if (!is.null(attr(path, 'status')) && attr(path, 'status')) {
    # user canceled
    return(NA)
  }

  #Error: Gtk-Message: GtkDialog mapped without a transient parent
  if(length(path) > 1){
    path = path[(length(path))]
  }

  return(path)
}

#' @name choose.dir.windows
#'
#' @title The windows version of the choose folder
#'
#' @seealso \code{\link{choose.dir}}
#'
#' @return
#' A length one character vector, character NA if 'Cancel' was selected.
#'
choose.dir.windows <- function(default = NA, caption = NA, useNew = TRUE) {
  if(useNew){
    ## uses a powershell script rather than the bat version, gives a nicer interface
    ## and allows setting of the default directory and the caption
    whereisutils <- system.file("utils", 'newFolderDialog.ps1', package = "shinyDirectoryInput")
    command = 'powershell'
    args = paste('-NoProfile -ExecutionPolicy Bypass -File',normalizePath(whereisutils))
    if (!is.null(default) && !is.na(default) && nzchar(default)) {
      args = paste(args, sprintf('-default "%s"', normalizePath(default)))
    }

    if (!is.null(caption) && !is.na(caption) && nzchar(caption)) {
      args = paste(args, sprintf('-caption "%s"', caption))
    }

    suppressWarnings({
      path = system2(command, args = args, stdout = TRUE)
    })
  } else {
    whereisutils <- system.file("utils", 'choose_dir.bat', package = "shinyDirectoryInput")
    command = normalizePath(whereisutils)
    args = if (is.na(caption)) '' else sprintf('"%s"', caption)
    suppressWarnings({
      path = system2(command, args = args, stdout = TRUE)
    })
  }
  if (path == 'NONE') path = NA
  return(path)
}


#' @name directoryInput
#'
#' @title Directory Selection Control
#'
#' @param inputId The \code{input} slot that will be used to access the value
#' @param label Display label for the control, or NULL for no label
#' @param value Initial value.  Paths are expanded via \code{\link{path.expand}}.
#'
#' @details
#' This widget relies on \code{\link{choose.dir}} to present an interactive
#' dialog to users for selecting a directory on the local filesystem.  Therefore,
#' this widget is intended for shiny apps that are run locally - i.e. on the
#' same system that files/directories are to be accessed - and not from hosted
#' applications (e.g. from shinyapps.io).
#'
#' @return
#' A directory input control that can be added to a UI definition.
#'
#' @seealso
#' \code{\link{updateDirectoryInput}}, \code{\link{readDirectoryInput}}, \code{\link{choose.dir}}
#'
directoryInput = function(inputId, label, value = NULL) {
  if (!is.null(value) && !is.na(value)) {
    value = path.expand(value)
  }
  version <- as.character(packageVersion("shinyDirectoryInput")[[1]])
  dep <- htmltools::htmlDependency(
    name = "shinyDirectoryInput-assets", version = version,
    package = "shinyDirectoryInput",
    src = "assets",
    script = "js/directory_input_binding.js"
  )
  tagList(
    shiny::div(
      class = 'form-group directory-input-container',
      shiny:::`%AND%`(label, tags$label(label)),
      shiny::div(
        shiny::span(
          class = 'col-xs-9 col-md-11',
          style = 'padding-left: 0; padding-right: 5px;',
          shiny::div(
            class = 'input-group shiny-input-container',
            style = 'width:100%;',
            div(class = 'input-group-addon', icon('folder-o')),
            tags$input(
              id = sprintf('%s__chosen_dir', inputId),
              value = value,
              type = 'text',
              class = 'form-control directory-input-chosen-dir',
              readonly = 'readonly'
            )
          )
        ),
        shiny::span(
          class = 'shiny-input-container',
          tags$button(
            id = inputId,
            class = 'btn btn-default directory-input',
            '...'
          )
        )
      )
    ),
    dep
  )

}


#' @name updateDirectoryInput
#'
#' @title Change the value of a directoryInput on the client
#'
#' @param session The \code{session} object passed to function given to \code{shinyServer}.
#' @param inputId The id of the input object.
#' @param value A directory path to set
#' @param ... Additional arguments passed to \code{\link{choose.dir}}.  Only used
#'    if \code{value} is \code{NULL}.
#'
#' @details
#' Sends a message to the client, telling it to change the value of the input
#' object.  For \code{directoryInput} objects, this changes the value displayed
#' in the text-field and triggers a client-side change event.  A directory
#' selection dialog is not displayed.
#'
updateDirectoryInput = function(session, inputId, value = NULL, ...) {
  if (is.null(value)) {
    value = choose.dir(...)
  }
  session$sendInputMessage(inputId, list(chosen_dir = value))
}

#' @name readDirectoryInput
#'
#' @title Read the value of a directoryInput
#'
#' @param session The \code{session} object passed to function given to \code{shinyServer}.
#' @param inputId The id of the input object
#'
#' @details
#' Reads the value of the text field associated with a \code{directoryInput}
#' object that stores the user selected directory path.
#'
readDirectoryInput = function(session, inputId) {
  session$input[[sprintf('%s__chosen_dir', inputId)]]
}
