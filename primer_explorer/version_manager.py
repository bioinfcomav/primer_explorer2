#!/usr/bin/env python3
import os
import subprocess


class VersionManager:
    version_fpath = os.path.join(os.path.dirname(os.path.realpath(__file__)), '_version.py')
    path = os.path.dirname(os.path.realpath(__file__))

    def _run_git_command(self, args):
        prefix = ['git', '-C', self.path]
        try:
            return subprocess.check_output(prefix + args).decode().strip()
        except subprocess.CalledProcessError:
            return None

    def get_git_version(self):
        git_version = self._run_git_command(['describe', '--tags', '--dirty=.dirty'])
        return git_version

    def get_file_version(self):
        with open(self.version_fpath, 'r') as fhand:
            version_in_file = fhand.readline().strip().split('=')[1]
            return version_in_file.strip()

    @property
    def version(self):
        git_version = self.get_git_version()
        file_version = self.get_git_version()
        if git_version != file_version:
            self.version = git_version
        return git_version

    @version.setter
    def version(self, new_version):
        with open(self.version_fpath, 'w') as fhand:
            fhand.write('version = "{}"\n'.format(new_version))
            fhand.flush()

    def update_version(self, pre_commit=True):
        git_version = self.get_git_version()
        items = git_version.split('-')
        if pre_commit:
            git_version = items[0]
            try:
                commit_num = int(items[1]) + 1
                git_version += '.dev{}'.format(commit_num)
            except IndexError:
                pass
        else:
            git_version = git_version.replace('-', '.dev', 1).replace('-', '+')[1:]
        self.version = git_version.replace('.dirty', '')


if __name__ == '__main__':
    import sys
    if sys.argv[1] == 'update_version':
        manager = VersionManager()
        manager.update_version(pre_commit=True)
