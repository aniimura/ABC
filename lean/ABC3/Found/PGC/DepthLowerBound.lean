import ABC3.Found.PGC.ExitLossDepth

/-!
# [pGC] ★②′/⑤ の射程は「無料の段」では変わらない —— `c d` は深さ `d` の**全部**を覆う

## ★どれを選んだか、なぜか（持ち場が選択を任せた点）

★**(D)（穴の再測定）を選んだ。** 4 候補の費用を先に測った:

| 候補 | 測った費用 | 判定 |
|---|---|---|
| (A) `a` の存在（不分岐 ⇒ 跡が全射） | `IsDiscreteValuationRing.TFAE` からの取り出し ＋ `IsIntegrallyClosed` / `Module.IsTorsionFree` / `Module.Finite`。★`grep -rn "IsIntegrallyClosed" lean/ABC3/Found/PGC/` は整数環については 0 件（在るのは `LubinTateCompletionDegree.lean:138` の別文脈） | ★重い |
| (B) `wildDepth` の接続 | `WildDepthDescent.lean:667 exists_natDegree_minpoly_descent_div` は★**`1/p` の平均に固定**されており（`:149` の表が「`p` で割った版」と書いている）、一般の `a` で書き直すには 1 ファイルぶん要る | ★重い |
| (C) `hform` | 新しい降下の構成そのもの | ★重い |
| ★(D) 穴の再測定 | 測定器 1 本（`AxWildDescent` を展開して `ε` で割るだけ） | ★**安い** |

## ★★★測定の結論 —— ②′/⑤ の射程は**変わらない**

前波の `UnramifiedStepFixed.exists_fixed_step_free` は
★「不分岐の段は `c = g = 1`（無料）」を示した。★**では②′（損失が定数 `> 1` で下から
抑えられれば積は非有界）や⑤（出口の損失は `axDecay p 1` 以上）は弱まるのか。**

★**弱まらない。** 理由は `AxWildDescent K c` の量化の順序である
（`AxTowerDecay.lean:473`）: `c` は `x` より**外**にあるので、
★`c d` は深さ `d` の**すべての** `x` を覆わねばならない。
⇒ ★**「ある `x` は無料で降りられる」ことは `c d` を下げない。**

これを測定器にしたのが §1 `le_c_of_lower_bound`:

★深さ `d ≠ 0` の `x` が 1 つでも「深さが下がるどの `x′` にも `D·ε` 以上離れている」なら
★**`D ≤ c d`**。

§2 `forall_le_c_of_witness` / `not_forall_prod_of_witness`:
★各深さ `d ≥ 1` にそういう `x`（証人）が在れば `∀ d ≥ 1, C ≤ c d` となり、
②′（`ExitLossDepth.not_forall_prod_Icc_le_of_le`）で積は非有界になる。

§3 `escape_needs_no_witness`（対偶）:
★★**②′ から逃げるには、ある深さ `d ≥ 1` で「証人が 1 つも無い」——
すなわち★深さ `d` の `x` が**全部**安く降りられる——ことが必要である。**

## ★まだ測っていないこと（正直に）

★**証人が各深さに在るか**は測っていない。
★見込みとしては `ℚ_p(p^{1/p^d})` のような完全分岐の `x` が証人になりそうだが、
★**それは形式化していないし、下界 `D·ε` の値も計算していない**。
⇒ ★②′/⑤ が「本当に効く」ことの最後の 1 点はここである。

## ★穴の現状（本波で射程の確認）

| 穴 | 本波の測定 |
|---|---|
| ①不分岐 | (b) は前波で閉じた。(a) が残る |
| ②′ | ★**射程は変わらない**（`c d` が全部を覆うため）。★ただし「証人の存在」が未測定 |
| ③幾何減衰 | 変わらず（跳びの変位の比の話で、`c` の量化とは無関係） |
| ⑤出口の限界 | ★**射程は変わらない**（②′と同じ理由） |
| ④ | 消えたまま |

★★**②は本日 2 度射程が動いた**（`c` → `g`、そして本波で「動かない」の確認）。
★動いていないことを確かめるのも測定である。

## ★①不分岐側との関係

★本波は①と**独立**である。測定器は分岐に一切触れていない
（`AxWildDescent` の定義を展開して `ε` で割るだけ）。

## 逸脱の記録（CLAUDE.md「逸脱」）

1. `.src` は上流と同じ pGC 物理 p.6 Corollary 3.1。
2. 既存の `Found/PGC/*.lean` は **1 行も書き換えていない**（import のみ）。
3. ★本ファイルは抽象核に切り出していない。理由: 中身が
   `AxWildDescent` の定義の展開そのもので、切り出すと核が空になる。
4. ★証人の存在を仮説（`hwit`）のまま残したことを上に明記した。
-/

namespace ABC3.Found.PGC

namespace DepthLowerBound

open ABC3.Skeleton.PGC

variable {p : ℕ} [Fact p.Prime]

/-! ## §1 測定器 —— 下界が 1 つあれば `c d` の下界になる -/

section Measure

/-- `wildDepth K x ≠ 0` なら `p ∣ deg minpoly x`。 -/
theorem dvd_natDegree_of_wildDepth_ne_zero (K : PAdicLocalField p) {x : K.closure}
    (h : wildDepth K x ≠ 0) : p ∣ (minpoly K.carrier x).natDegree := by
  by_contra hd
  exact h ((wildDepth_eq_zero_iff K x).mpr hd)

/-- ★★★**測定器** —— 深さ `d ≠ 0` の `x` が 1 つでも
「深さが下がるどの `x′` にも `D·ε` 以上離れている」なら `D ≤ c d`。

★`AxWildDescent K c`（`AxTowerDecay.lean:473`）は `c` を `x` より**外**で量化するので、
★`c d` は深さ `d` の**すべての** `x` を覆わねばならない。
⇒ ★「ある `x` は無料で降りられる」ことは `c d` を下げない。 -/
theorem le_c_of_lower_bound (K : PAdicLocalField p) {c : ℕ → ℝ}
    (h : AxWildDescent K c) {ε D : ℝ} (hε : 0 < ε) {x : K.closure}
    (hdep : wildDepth K x ≠ 0) (hx : ∀ σ : K.absGal, ‖σ • x - x‖ ≤ ε)
    (hlow : ∀ x' : K.closure, wildDepth K x' < wildDepth K x → D * ε ≤ ‖x - x'‖) :
    D ≤ c (wildDepth K x) := by
  obtain ⟨x', hlt, hd, _⟩ :=
    h ε (le_of_lt hε) x hx (dvd_natDegree_of_wildDepth_ne_zero K hdep)
  have h1 : D * ε ≤ c (wildDepth K x) * ε := le_trans (hlow x' hlt) hd
  exact le_of_mul_le_mul_right (by linarith) hε

end Measure

/-! ## §2 各深さに証人が在れば②′がそのまま効く -/

section Scope

/-- ★各深さ `d ≥ 1` に証人が在れば `∀ d ≥ 1, C ≤ c d`。 -/
theorem forall_le_c_of_witness (K : PAdicLocalField p) {c : ℕ → ℝ} {C : ℝ}
    (h : AxWildDescent K c)
    (hwit : ∀ d : ℕ, 1 ≤ d → ∃ (ε : ℝ) (x : K.closure), 0 < ε ∧ wildDepth K x = d ∧
      (∀ σ : K.absGal, ‖σ • x - x‖ ≤ ε) ∧
      (∀ x' : K.closure, wildDepth K x' < d → C * ε ≤ ‖x - x'‖)) :
    ∀ d : ℕ, 1 ≤ d → C ≤ c d := by
  intro d hd
  obtain ⟨ε, x, hε, hdx, hx, hlow⟩ := hwit d hd
  have hdep : wildDepth K x ≠ 0 := by rw [hdx]; omega
  have := le_c_of_lower_bound K h hε hdep hx (by rw [hdx]; exact hlow)
  rwa [hdx] at this

/-- ★★②′（`ExitLossDepth.not_forall_prod_Icc_le_of_le`）がそのまま効く形。
★前波の「不分岐の段は無料」は②′を弱めない。 -/
theorem not_forall_prod_of_witness (K : PAdicLocalField p) {c : ℕ → ℝ} {C B : ℝ}
    (hC : 1 < C) (h : AxWildDescent K c)
    (hwit : ∀ d : ℕ, 1 ≤ d → ∃ (ε : ℝ) (x : K.closure), 0 < ε ∧ wildDepth K x = d ∧
      (∀ σ : K.absGal, ‖σ • x - x‖ ≤ ε) ∧
      (∀ x' : K.closure, wildDepth K x' < d → C * ε ≤ ‖x - x'‖)) :
    ¬ (∀ n : ℕ, ∏ d ∈ Finset.Icc 1 n, c d ≤ B) :=
  ExitLossDepth.not_forall_prod_Icc_le_of_le hC (forall_le_c_of_witness K h hwit)

end Scope

/-! ## §3 ★★逃げるには「ある深さで証人が 1 つも無い」ことが要る -/

section Escape

/-- ★★★**逃げ道の必要条件** —— `c d < C` なら、深さ `d` には証人が 1 つも無い
（すなわち深さ `d` の `x` が**全部**安く降りられる）。★`d = 0` は自明な場合。 -/
theorem escape_needs_no_witness (K : PAdicLocalField p) {c : ℕ → ℝ} {C : ℝ} {d : ℕ}
    (h : AxWildDescent K c) (hlt : c d < C) :
    ¬ (∃ (ε : ℝ) (x : K.closure), 0 < ε ∧ wildDepth K x = d ∧
        (∀ σ : K.absGal, ‖σ • x - x‖ ≤ ε) ∧
        (∀ x' : K.closure, wildDepth K x' < d → C * ε ≤ ‖x - x'‖)) ∨ d = 0 := by
  rcases Nat.eq_zero_or_pos d with hd0 | hd0
  · exact Or.inr hd0
  · refine Or.inl ?_
    rintro ⟨ε, x, hε, hdx, hx, hlow⟩
    have hdep : wildDepth K x ≠ 0 := by rw [hdx]; omega
    have hle := le_c_of_lower_bound K h hε hdep hx (by rw [hdx]; exact hlow)
    rw [hdx] at hle
    linarith

end Escape

/-! ## §4 `.src` と 使っている公理の一覧 -/

section Src

def le_c_of_lower_bound.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

def forall_le_c_of_witness.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

def not_forall_prod_of_witness.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

def escape_needs_no_witness.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

end Src

#print axioms dvd_natDegree_of_wildDepth_ne_zero
#print axioms le_c_of_lower_bound
#print axioms forall_le_c_of_witness
#print axioms not_forall_prod_of_witness
#print axioms escape_needs_no_witness

end DepthLowerBound

end ABC3.Found.PGC
