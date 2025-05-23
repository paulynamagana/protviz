---
layout: default
title: Examples
nav_order: 4
description: "Understand the key examples."
---
# Advanced usage & examples
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


Beyond visualising single data types, the real power of Protviz comes from its ability to integrate and display information from multiple sources simultaneously. By layering various tracks, you can build a comprehensive picture of your protein's features, structural characteristics, and predicted properties.
{: .fs-6 .fw-300 }

This section provides practical examples of how to combine different data retrieval clients and track types to create richer visualisations. Each example will focus on a specific combination of data, offering a complete script that you can adapt for your own proteins of interest.
