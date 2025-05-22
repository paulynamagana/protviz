---
layout: default
title: Core Components
nav_order: 3
description: "Understand the key building blocks of Protviz: Data Retrieval Clients and Visualisation Tracks."
---

# Core Components: Data Retrieval Clients and Visualisation Tracks
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
Protviz is built around two key types of components: **Data Retrieval Clients** that fetch information from various bioinformatics databases, and **Visualisation Tracks** that display this information. Understanding these components will help you leverage the full power of the package.
{: .fs-6 .fw-300 }
