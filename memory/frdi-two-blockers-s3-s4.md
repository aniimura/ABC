---
name: frdi-two-blockers-s3-s4
description: [FrdI] §3 9/9 は Gap_1_14_iii(反例つき原典の穴)で塞がれ、§4 10/10 は pairing-vanishes(5 件の下流)で塞がれている。Prop 4.1 は 2026-08-18 に閉じた(§4 4/10)ので、到達可能な残りは Remark 4.5.1 の 1 件。
metadata:
  type: project
---

★2026-08-18 の測定。`Proposition 3.2` を閉じた(31/54、§3 8/9)あとの現在地。

## ★★§3 9/9 —— `Theorem 3.4` が**原典の穴**を消費する

- `Theorem 3.4, (i)` は**完全に実装済み**(`thm_3_4_i.src`、`Thm34.lean` ＋ `Thm34Rigid.lean`)。
- 残る **(ii)** は `Proposition 1.14, (ii)(iii)` の**両向き**を要求する
  (Ψ が pre-step を保つことを `BoundedFSMIFactor` という圏論的条件の同値で移すため)。
- ★`Prop 1.14, (iii)` の `⟸`(`prop_1_14_iii_mpr`)は **`Definition 1.3` から出ない** ——
  `Check/FrdI/TwistedFrobenioid.lean` に**機械検証済みの反例**(捻れ積
  `𝔽_ℕ ⋉ ∏ ℤ/n`)があり、`Gap/FrdI/Section1.lean` の `Gap_1_14_iii` に
  ③ `sourceGap` として登録済み。埋めるには原典の語彙 2 語
  (**`unit-trivial` 型** ← mono、**`Φ` が `perfect`** ← fiberwise-surjectivity)が要る。
- ★**迂回路は 3 本とも塞がった**(2026-08-18 に検討): `Prop 1.7, (iv)` 経由・
  準逆関手経由・既約 pre-step で鎖を作る。鎖を**非有界**にできるのは
  **大きな素数次数の prime-Frobenius 射**だけで、そこがまさに穴である。
- ★★したがって §3 9/9 は**原典に忠実なままでは届かない**。

## ★★§4 10/10 —— 5 件が `pairing-vanishes` の下流

`node tools/frdi-blocked.mjs`: §4 で**チェーンに依らず届くのは 5/10**。
`Thm 4.2` / `Prop 4.4` / `Prop 4.8` / `Cor 4.10` / `Cor 4.11` は
`Prop 4.4, (ii)` の `otriBase` = `𝒪^▷(X^birat)` の可換性待ち。

★★**`otriBase` は可換性と同値である**(緩められない)—— `birat_otriBase_of_comm` の
`h : φ ≫ β = α ≫ φ` は `β = φ⁻¹αφ` を意味し、`φ` は `𝒞^birat` で同型なので
`β` は `𝒪^▷(B)` 全体を走る。よって `otriBase ⟺ w = φ⁻¹φ' が全体と可換`。

★2026-08-18 の前進: スライス `(𝒞^coa-pre)_A` は**前順序**(`coaPreOverFunctor` が忠実で
行き先 `Order` の hom が subsingleton)なので、**持ち上げは一意**。これで Ore の脚の
一致 `α ≫ p' = β ≫ p` は出る。★残るのは `L(k≫f) ≫ g = L(k≫g) ≫ f`
(`L` = `c` に沿った一意の持ち上げ)——これが `pairing-vanishes` の正体。

## ★到達可能な残りは 1 件(§4 4/10 → 5/10)

- ★★`Proposition 4.1` は **2026-08-18 に (i)〜(v) 全部閉じた**。
  ★鍵だったのは、**(iv) と (v) の単系層が同一の 1 本**だと気づいたこと
  (`isPrimaryElt_iff_exists_prime_cond`)。原文の
  「reversing the direction of the arrows」は、圏の側で
  **後置(スライス)を前置(コスライス)に取り替えるだけ**になった。
  ★素点ごとの比較の実体は `mle_of_restrict`(`Def24SuppElt.lean`)。
- `Remark 4.5.1` —— standard 型は済(`istr_standardType`)。
  rationally standard 型の 3 条(birat Frobenius-normalized・rational・
  `(𝒞^un-tr)^birat` の Frobenius-compact 対象)が残り、いずれも
  `𝒞^istr` の birationalization / unit-trivialization を経由する。

**How to apply:** §3/§4 の残りを見積もるときは、まずこの 2 つの障害を確認する。
数字は `node tools/frdi-progress.mjs` と `node tools/frdi-blocked.mjs`。
関連: [[frdi-three-chains]] [[pf-frobenius-absorbed-by-transition]]
[[no-wall-decompose-instead]]
