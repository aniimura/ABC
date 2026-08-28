/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.PrimesOfSize
import ABC3.Found.GenEll.Elementary
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★★★★★★★`Corollary 4.3` の算術の鎖（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.22–p.23。

原文 (GenEll p.22):
> Corollary 4.3.

## ★★★★★★★★★★★★★★★★★★★これは何か

原文 p.22–p.23 の `Corollary 4.3` の証明は、`Theorem 3.8` を当てたあと
**素数の存在（`Lemma 4.1`）と不等式の帳簿**でできている。
★★★**本ファイルはその帳簿を取る。**

### ★原文の帳簿（p.22 末〜p.23）

    `x_{S∘} ≤ x_S + (1 + 3/2)·23040d·deg∞([E_L])`            （`Lemma 4.2`）
           `≤ x_S + 3·12·23040d·ht^Falt([E_L]) + d·C″`        （`Prop 3.4`、`1+ϵ = 6/5`）

    `x_{S•} ≤ x_S + 3·12·23040d·ht^Falt + 3d·log-diff_Mell + d·C″`

そして `Lemma 4.1` を `M = 1`、`1+6ϵ = 2`、`h = 23040·100d·(ht^Falt + C′·d^ϵ)`、
`A = S∘` または `S•` で当てると `h ≤ l ≤ 2·x_A + 8h` が出る。★原文の

> [and applying the estimate 2 · 3 · 12 + 8 · 100 ≤100 + 800 = 900]

がこの帳簿の要点である（`2·36 + 800 = 872 ≤ 900`）。

## ★★★何を受けているか（測定、2026-08-29）

| 入力 | 状態 |
|---|---|
| `Lemma 4.1`（大きさの指定された素数の存在） | ★`Found/GenEll/PrimesOfSize.lean` にある |
| `Lemma 4.2`（初等評価） | ★`Found/GenEll/Elementary.lean` にある |
| 素数定理 `θ(x)/x → 1` | ★`Found/GenEll/PrimeNumberTheorem.lean` にある |
| `Proposition 3.4`（`deg∞ ≲ 12(1+ϵ)ht^Falt`） | △ 第 1・第 2 の合成は `§9-980`、第 3 が残る |
| `Theorem 3.8`（Galois 像が `SL₂` を含む） | ☆残る（Serre の開像定理） |
| 「`l` が `L/ℚ` で不分岐なら `SL₂ ⊆ 像 ⟺ 全射`」 | ☆残る（`ℚ(ζ_{l^∞})/ℚ` が `l` で完全分岐、ゆえに `L/ℚ` と線型無関連） |
| ★不等式の帳簿 | ★★**本ファイル** |

★★★★`Corollary 4.3` に残るのは**`Theorem 3.8`（Serre の開像定理）と `Prop 3.4` の第 3 の `≲`**である。
★`.src` は条つき——指標には数えない。
-/

namespace ABC3.Found.GenEll

/-! ## ★★★★原文が名指しする数値の評価 -/

/-- ★★★★**原文の `2 · 3 · 12 + 8 · 100 ≤ 100 + 800 = 900`**。

原文 (GenEll p.22):
> Corollary 4.3.

★`2·36 + 800 = 872 ≤ 900`——`Lemma 4.1` の `2·x_A + 8h` に
`x_A` と `h` の評価を入れたときの係数である。 -/
theorem cor_4_3_coeff : 2 * (3 * 12) + 8 * 100 ≤ 100 + 800 := by norm_num

/-! ## ★★★★★★★★★★★★★★`l∘` の評価 -/

/-- ★★★★★★★★★★★★★★**`Corollary 4.3, (c)` の `l∘` の評価**。

原文 (GenEll p.22):
> Corollary 4.3.

    `l∘ ≤ 23040 · 900d · ht^Falt([E_L]) + 2x_S + C · d^{1+ϵ}`

★仮定 `hxSo`（`Lemma 4.2`）・`hdeg`（`Prop 3.4`、`1+ϵ = 6/5`）・
`hl`（`Lemma 4.1` を `M = 1`、`1+6ϵ = 2`、`h = 23040·100d·(ht^Falt + C′d^ϵ)` で当てた形）から、
**係数 `2·3·12 + 8·100 = 872 ≤ 900`** で結論が出る。

★★`B ≤ ht^Falt`（`ht^Falt` が下に有界であること、`Prop 3.4`）を使って
`872 ≤ 900` の差を `|B|·d^{1+ϵ}` に吸収している——★原文が `C` に押し込んでいる分である。 -/
theorem cor_4_3_arith (d xS xSo dgInf htF B C' C'' eps : ℝ)
    (hd : 1 ≤ d) (heps : 0 < eps) (hC'' : 0 ≤ C'')
    (hB : B ≤ htF)
    (hxSo : xSo ≤ xS + (5 / 2) * (23040 * d) * dgInf)
    (hdeg : dgInf ≤ (12 * (6 / 5)) * htF + C'')
    (l : ℝ) (hl : l ≤ 2 * xSo + 8 * (23040 * 100 * d * (htF + C' * d ^ eps))) :
    l ≤ 23040 * 900 * d * htF + 2 * xS
        + (23040 * (5 * C'' + 800 * C' + 28 * |B|)) * d ^ (1 + eps) := by
  have hdpos : (0:ℝ) < d := by linarith
  have hsplit : d ^ (1 + eps) = d * d ^ eps := by
    rw [Real.rpow_add hdpos, Real.rpow_one]
  have hdle : d ≤ d ^ (1 + eps) := by
    calc d = d ^ (1:ℝ) := (Real.rpow_one d).symm
      _ ≤ d ^ (1 + eps) := Real.rpow_le_rpow_of_exponent_le hd (by linarith)
  have hdegd : (23040 * d) * dgInf
      ≤ (23040 * d) * ((12 * (6 / 5)) * htF + C'') :=
    mul_le_mul_of_nonneg_left hdeg (by positivity)
  have hBd : 28 * (23040 * d) * B ≤ 28 * (23040 * d) * htF :=
    mul_le_mul_of_nonneg_left hB (by positivity)
  have hBabs : -(28 * 23040 * |B|) * d ^ (1 + eps) ≤ 28 * (23040 * d) * B := by
    have h1 : -|B| ≤ B := neg_abs_le B
    nlinarith [hdle, abs_nonneg B]
  have hC''d : (5 * 23040) * d * C'' ≤ (5 * 23040) * C'' * d ^ (1 + eps) := by
    nlinarith [hdle, hC'']
  have hC'd : 800 * 23040 * (d * d ^ eps) * C' = 800 * 23040 * C' * d ^ (1 + eps) := by
    rw [hsplit]; ring
  nlinarith [hl, hxSo, hdegd, hBd, hBabs, hC''d, hC'd]

/-! ## ★★★★★★★★★★★★★★`l•` の評価（`log-diff` つき） -/

/-- ★★★★★★★★★★★★★★**`Corollary 4.3, (c)` の `l•` の評価**。

原文 (GenEll p.22):
> Corollary 4.3.

    `l• ≤ 23040 · 900d · ht^Falt + 6d · log-diff_Mell + 2x_S + C · d^{1+ϵ}`

★`S•` は `S∘` に「`L/ℚ` で分岐する素数」と「分岐指数を割る素数」を足したものであり、
原文はその寄与を `3d·log-diff_Mell` で押さえる
（『the primes appearing in the arithmetic divisor that gives rise to log-diff_Mell
appear with multiplicity ≥ one less than the ramification indices of L/ℚ』）。
★★`Lemma 4.1` の `2·x_A` で **2 倍**になり `6d·log-diff` になる。 -/
theorem cor_4_3_arith_bullet (d xS xSb ldiff htF B C' C'' eps : ℝ)
    (hd : 1 ≤ d) (heps : 0 < eps) (hC'' : 0 ≤ C'')
    (hB : B ≤ htF)
    (hxSb : xSb ≤ xS + (3 * 12 * 23040) * d * htF + 3 * d * ldiff + d * C'')
    (l : ℝ) (hl : l ≤ 2 * xSb + 8 * (23040 * 100 * d * (htF + C' * d ^ eps))) :
    l ≤ 23040 * 900 * d * htF + 6 * d * ldiff + 2 * xS
        + (23040 * (2 * C'' + 800 * C' + 28 * |B|)) * d ^ (1 + eps) := by
  have hdpos : (0:ℝ) < d := by linarith
  have hsplit : d ^ (1 + eps) = d * d ^ eps := by
    rw [Real.rpow_add hdpos, Real.rpow_one]
  have hdle : d ≤ d ^ (1 + eps) := by
    calc d = d ^ (1:ℝ) := (Real.rpow_one d).symm
      _ ≤ d ^ (1 + eps) := Real.rpow_le_rpow_of_exponent_le hd (by linarith)
  have hBd : 28 * (23040 * d) * B ≤ 28 * (23040 * d) * htF :=
    mul_le_mul_of_nonneg_left hB (by positivity)
  have hBabs : -(28 * 23040 * |B|) * d ^ (1 + eps) ≤ 28 * (23040 * d) * B := by
    have h1 : -|B| ≤ B := neg_abs_le B
    nlinarith [hdle, abs_nonneg B]
  have hC''d : 2 * d * C'' ≤ 2 * C'' * d ^ (1 + eps) := by
    nlinarith [hdle, hC'']
  have hC'd : 800 * 23040 * (d * d ^ eps) * C' = 800 * 23040 * C' * d ^ (1 + eps) := by
    rw [hsplit]; ring
  nlinarith [hl, hxSb, hBd, hBabs, hC''d, hC'd]

/-! ## ★出典の紐付け(`.src`)——★**条つきである。指標には数えない** -/

def cor_4_3_coeff.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 22,
    item := "Corollary 4.3(原文の 2·3·12 + 8·100 ≤ 100 + 800 = 900)",
    sectionId := "genell-cor-4-3" }

def cor_4_3_arith.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 22,
    item := "Corollary 4.3, (c)(l∘ の評価——不等式の帳簿)",
    sectionId := "genell-cor-4-3" }

def cor_4_3_arith_bullet.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 22,
    item := "Corollary 4.3, (c)(l• の評価——log-diff つき)",
    sectionId := "genell-cor-4-3" }

def cor_4_3_arith.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "lemma_4_1(大きさの指定された素数の存在)"
      (.inProject "ABC3" "ABC3.Found.GenEll.lemma_4_1") 4,
    .citation "[ABC3]" "lemma_4_2(初等評価)"
      (.inProject "ABC3" "ABC3.Found.GenEll.lemma_4_2") 2,
    .otherPaper "[GenEll]"
      ("Proposition 3.4(deg∞ ≲ 12(1+ϵ)·ht^Falt、ϵ を 1/5 に取る)" ++
       "——★第 1・第 2 の合成は §9-980 で取れた。第 3 の ≲ が残る") 8,
    .otherPaper "[GenEll]"
      "Theorem 3.8(Galois 像が SL₂ を含む)——★Serre の開像定理。**残る**" 11,
    .folklore
      ("『l が L/ℚ で不分岐なら SL₂ ⊆ 像 ⟺ 全射』" ++
       "——ℚ(ζ_{l^∞})/ℚ が l で完全分岐、ゆえに L/ℚ と線型無関連。**残る**") 5,
    .implicitStep
      ("★★★★★★測定(2026-08-29): Corollary 4.3 の証明は Theorem 3.8 を当てたあと" ++
       "**素数の存在(Lemma 4.1)と不等式の帳簿**でできている。本ファイルはその帳簿を取る。" ++
       "★原文が名指しする係数 2·3·12 + 8·100 ≤ 100 + 800 = 900 がその要点であり、" ++
       "872 ≤ 900 の差は ht^Falt が下に有界であること(B ≤ ht^Falt)を使って" ++
       "|B|·d^{1+ϵ} に吸収している——原文が C に押し込んでいる分である") 7,
    .implicitStep
      ("★★★Corollary 4.3 に残るのは **Theorem 3.8(Serre の開像定理)と " ++
       "Proposition 3.4 の第 3 の ≲** である。" ++
       "★Lemma 4.1・Lemma 4.2・素数定理は本プロジェクトにある" ++
       "(Found/GenEll/PrimesOfSize.lean／Elementary.lean／PrimeNumberTheorem.lean)") 7 ]

end ABC3.Found.GenEll
