document.querySelectorAll('form.pick').forEach(form => {
  const order = JSON.parse(form.dataset.order);
  const root = form.dataset.root;
  const boxes = [...form.querySelectorAll('input[type=checkbox]')];
  const btn = form.querySelector('button');
  const hint = form.querySelector('.hint');
  let picked = boxes.filter(b => b.checked).map(b => b.value);
  function update() {
    btn.disabled = picked.length === 0;
    hint.textContent = picked.length === 0 ? 'Pick one microbiome to explore, or two to compare.'
      : picked.length === 1 ? 'Add a second microbiome to compare, or show this one alone.' : '';
  }
  boxes.forEach(b => b.addEventListener('change', () => {
    if (b.checked) {
      picked.push(b.value);
      if (picked.length > 2) {   // keep the two most recent choices
        const drop = picked.shift();
        boxes.find(x => x.value === drop).checked = false;
      }
    } else {
      picked = picked.filter(v => v !== b.value);
    }
    update();
  }));
  form.addEventListener('submit', e => {
    e.preventDefault();
    const id = [...picked].sort((a, b) => order.indexOf(a) - order.indexOf(b)).join('--');
    location.href = root + 'views/' + id + '/index.html';
  });
  update();
});
