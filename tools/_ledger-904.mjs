// 2026-08-31（第 899-904）—— Lemma 3.5 の残りを台帳に書き直す。
// ★CLAUDE.md の規約により、台帳の更新は node -e ではなく .mjs で書く。
import { readFileSync, writeFileSync } from 'node:fs';

const P = 'D:/Math_ABC3/ResearchPaper/mathlib-gap.json';
const d = JSON.parse(readFileSync(P, 'utf8'));

d.lemma35Remaining20260831 = {
  measuredAt: '2026-08-31',
  block: '第 899-904',
  metric: 'GenEll 必要分 18 / 24（§1 9/9・§2 0/1・§3 6/9・§4 3/5）——指標は動いていない',
  closedThisStretch: [
    '★Tate 一意化 Φ の存在（第 899 `dvrTatePhiAddEquiv`）——' +
      'tatePhiAddEquivAll の 6 つの仮説はすべて DVR の事実（局所環・極大は素・整閉・付値環）から作れる。',
    '★★★★Néron–Ogg–Shafarevich の義務（第 902）——' +
      'Lemma 3.5 の証明が使うのは `l·deg∞(E) ≤ deg∞(E′)` の片側だけであり、' +
      '良い素点では左辺が 0 なので `minDeltaExp_nonneg` だけで自動的に成り立つ。' +
      'したがって `jExp_velu_good`（同種な曲線は同じ還元型をもつ）は**要らない**。',
    '★局所の連鎖の終点（第 904 `minDeltaExp_eq_mul_of_veluMu`）——' +
      '数体の素点の完備化で `v_p(Δ_min(E′)) = l·v_p(Δ_min(E))` が出る。',
  ],
  remaining: [
    '☆(1) 分裂乗法還元: `SemistableAt p E` ＋ `jExp p E < 0` から ' +
      '`(E ⊗ Lv).HasSplitMultiplicativeReduction R` を作る。' +
      '★mathlib の HasSplitMultiplicativeReduction は剰余体で 2 次式が分裂することを要求するので、' +
      '非分裂の場合は不分岐 2 次拡大へ上げる段が要る。',
    '☆(2) 極小モデルの移行: `IsMinimal (primeSubring p) (C • E)` から ' +
      '`IsMinimal R (E ⊗ Lv)` を作る（極小判別式の付値は完備化で変わらない）。',
    '☆(3) H の像: `l` が局所高さと互いに素なら l-巡回部分群 H は ' +
      '`⟨Φ(ζ)⟩`（= μ_l）に対応する。★`lemma_3_2_i_tate_all`（証明済み）の対偶。',
  ],
  note:
    '★これ以外の数学はすべて証明済みである（第 704・853-904）。' +
    '外部文献への未証明の依存は無い（ht^Falt は第 704 で内製、NOS は第 902 で不要になった）。',
};

writeFileSync(P, JSON.stringify(d, null, 2) + '\n', 'utf8');
console.log('mathlib-gap.json: lemma35Remaining20260831 を書いた');
