---
title: Running Examples
nav_order: 5
---

# Running Examples
{: .no_toc }
<button class="btn js-toggle-dark-mode">Dark mode</button>

<script>
const toggleDarkMode = document.querySelector('.js-toggle-dark-mode');

jtd.addEvent(toggleDarkMode, 'click', function(){
  if (jtd.getTheme() === 'dark') {
    jtd.setTheme('light');
    toggleDarkMode.textContent = 'Dark mode';
  } else {
    jtd.setTheme('dark');
    toggleDarkMode.textContent = 'Light mode';
  }
});
</script>
---

The package includes several example scripts (e.g., example_pdbe.py, example_afdb.py) in the examples/ directory.

To run an example, navigate to the directory containing the scripts and execute it with Python:

```bash
python examples/example_afdb.py
```
