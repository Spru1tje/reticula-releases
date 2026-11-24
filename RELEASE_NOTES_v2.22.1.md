# Reticula 2.22.1

Hotfix release.

- **Agglomerate tagger on the TEM page**: an expander "Agglomerate tag: rods and clumps (first image)" lets you click
  objects and tag them agglomerate or rod. The tagger refits its model on all tags and recounts the picture.
  Segmentation is cached, so a click does not segment again. A fault in the tagger is shown in the expander only.
- Release fix: the 2.22.1 tag had no release notes file, so the release workflow stopped at the prepare step.

Known limits:
- The tagger has no real tags yet; its accuracy on your micrographs is unknown.
- Precise mode and rod measurement have no manual reference yet.

Full list: app/CHANGELOG.md.
