#For dealing with file paths properly
from pathlib import Path
#Type hints
from typing import *
#Enum support
from enum import Enum, auto


#"Path requirements" - enum used to describe the requisite state of a path being
#added to the ledger
class PathReq(Enum):
	#No special requirements
	NONE     = auto()
	#Path must point to existing file
	IS_FILE = auto()
	#Path must not point to existing file
	NO_FILE  = auto()
	#Path must point to existing directory
	IS_DIR  = auto()
	#Path must not point to existing directory
	NO_DIR   = auto()

def pathchk(path: Path, behavior: PathReq = PathReq.NONE) -> None:
	"""Determines whether a path fulfills a given requirement.
	
	Args:
		path: pathlib.Path
			The path to check.
		behavior: PathReq
			The requirement for the path. See PathReq class for description.
			Defaults to PathReq.None if out of range.
	
	Raises:
		FileNotFoundError
		FileExistsError
			Depending on the on-disk state of `path' and the requested behavior.
	"""
	
	match behavior:
		case PathReq.IS_FILE:
			if not path.is_file():
				raise FileNotFoundError(f"File {path} does not exist")
		
		case PathReq.NO_FILE:
			if path.is_file():
				raise FileExistsError(f"File {path} already exists")
		
		case PathReq.IS_DIR:
			if not path.is_dir():
				raise FileNotFoundError(f"Directory {path} does not exist")
		
		case PathReq.NO_DIR:
			if path.is_dir():
				raise FileExistsError(f"Directory {path} already exists")


#Class acting as a ledger for files and directories
class PathLedger:
	def __init__(
			self
		) -> None:
		"""Instantiates a PathLedger."""
		
		#Dictionary of items in ledger.
		#Keys are common names, values are actual paths.
		self.__dict__["d"]: dict[str:Path] = dict()
	
	
	def __setattr__(
			self,
			name:	str,
			params:	str|Path|tuple[str|Path]|tuple[str|Path,int]|list
		) -> None:
		"""Adds a new item to the ledger. Can be either a file, or a directory.
		
		Args:
			name: str
				The new item's common name.
			params:	str|Path|tuple[str|Path]|tuple[str|Path,int]|list
				Two parameters describing the creation of the new item.
				__setattr__() can only accept three total arguments, so these
				parameters have to be packed, either into a tuple or a list,
				unless only the first, required parameter is passed.
				
				The first/only value (str or pathlib.Path) is the path of the
				new item.
				
				The second item (int) is an optional requirement that the path
				should fulfill. If it does not fulfill the requirement, an
				exception will be raised
		
		Raises:
			FileNotFoundError
			FileExistsError
				If the path does not fulfill the passed requirement.
			TypeError
				If `params' was an inappropriate type.
		"""
		
		#Deencapsulate `params' as needed
		if isinstance(params, (str, Path)):
			#If a str was passed, encapsulate it to a pathlib.Path
			#Calling pathlib.Path() on a pathlib.Path() instance won't
			#re-encapsulate it
			path: Path = Path(params)
			#Default behavior value
			behavior: PathReq = PathReq.NONE
		
		elif isinstance(params, (tuple, list)):
			#If a tuple/list was passed, unpack it
			#Encapsulate the path to a proper pathlib.Path as needed
			path: Path = Path(params[0])
			behavior: PathReq = params[1] if len(params) > 1 else PathReq.NONE
		
		else:
			#Refuse to deal with improperly packed params
			raise TypeError(f"Bad parameters for new ledger item '{name}'")
		
		
		#Add new path to dictionary under the passed common name
		self.d[name] = path
		
		
		#Check if path fulfills requirement, if any
		pathchk(path, behavior)
	
	
	def __getattr__(
			self,
			name:	str
		) -> Path:
		"""Returns the path of an item in the ledger, based on its common name.
		
		Args:
			name: str
				Name of the item to return.
		
		Returns:
			pathlib.Path
				Path of the requested item.
		
		Raises:
			KeyError
				If the requested item is not in the ledger.
		"""
		if name in self:
			return self.d[name]
		else:
			raise KeyError(f"Item '{name}' not found in path ledger")
	
	
	def __contains__(
			self,
			name:	str
		) -> bool:
		"""Checks whether or not an item of the given name is in the ledger.
		Minimal shortcut over having to address the dictionary directly.
		
		Args:
			name: str
				Name of the item whose presence to check for.
		
		Returns:
			bool
				True if an item of the given name is in the ledger.
		"""
		return name in self.d
