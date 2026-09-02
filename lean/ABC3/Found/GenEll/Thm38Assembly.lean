/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Skeleton.GenEll.Section3
import ABC3.Found.GenEll.Lemma31
import ABC3.Meta.Claim

/-!
# 第 1445 ブロック —— **★★★★★★★★★★★★★★★★
`Theorem 3.8` の証明本体を `Found/` へ**

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19–p.20。

原文 (GenEll p.19):
> Theorem 3.8. (Full Special Linear Galois Actions) Let Q[bb][bar] be an algebraic

## ★★これは何か——第 1444 と同じ、場所を規約どおりに直す仕事である

`Theorem 3.8` の証明は `sorry` 0 であったが、本体が
`Skeleton/GenEll/GaloisImage.lean` に置かれたままであった。
☆規約は `Skeleton/` = 原典の主張を写す場所、`Found/` = `sorry` の無い実装であり、
`Lemma 4.1` / `Lemma 4.2` / `Corollary 4.3` / `Corollary 4.4`（第 1444）は
既にその形（`Skeleton` は写して `Found` へ委譲）になっている。

## ★★★★中身は 1 文字も変えていない

`#print axioms` は移動の前後どちらでも
`[propext, Classical.choice, Quot.sound]` である（`sorryAx` は無い）。

## ☆筋

`lemma_3_7` から `C₇`・`Exc` を取り、`E′ = torsionExt E` へ移って
`imageContainsSL2_of_torsionExt` を当てる。
★**`5 ≤ l` はどちらの枝からも出る**（第 1434、`Check/GenEll/VeluQuotOKNeedsL5.lean` が機械検証）:
枝 (a) は括弧の中が `≥ 1` なので `l ≥ 2304000`、
枝 (b) は `Nat.Coprime l 30` から `interval_cases` で出る。

## ★★欄は `Interface` で受けたままである

`EllModuliData` の欄（`imageContainsSL2_of_torsionExt` 等）は未構成であり、
これは「`Theorem 3.8` が絶対的に証明された」ことを意味しない。
☆欄を埋めるのは `Skeleton/GenEll/EllModuliWitness.lean` の仕事で、そこには
`VeluSemistable.lean` の `j(E′)` の整性と `GaloisLocal.lean` の Tate 加群の Galois 像が残っている。
-/

namespace ABC3.Found.GenEll

open ABC3.Meta ABC3.Interface.GenEll ABC3.Skeleton.GenEll

set_option maxHeartbeats 1600000 in
/-- **[GenEll] Theorem 3.8**(Full Special Linear Galois Actions)——
★**証明本体**（第 1445 で `Skeleton` から移した）。

原文 (GenEll p.19):
> Theorem 3.8. (Full Special Linear Galois Actions) Let Q[bb][bar] be an algebraic

☆statement の逐語の読み（`ϵ` を `.txt` から読んではならないこと、`23040` の由来）は
`Skeleton/GenEll/GaloisImage.lean` の `theorem_3_8` にある。 -/
theorem theorem_3_8 (D : EllModuliData)
    (KV : Set D.EllClass) (hKV : D.CompactlyBounded KV)
    (ε : ℝ) (hε : 0 < ε) :
    ∃ C : ℝ, 0 < C ∧ ∃ Exc : Set D.EllClass, D.GaloisFinite Exc ∧
      ∀ (E : D.Curve) (l : ℕ), Nat.Prime l → D.cls E ∉ Exc →
        ((23040 * 100 * (D.degOfDefinition E : ℝ)
              * (D.faltingsHeight (D.cls E) + C * (D.degOfDefinition E : ℝ) ^ ε) ≤ (l : ℝ)
            ∧ D.HasPotMultRed E)
          ∨ (D.cls E ∈ KV ∧ D.PrimeToLocalHeights E l ∧ Nat.Coprime l 30)) →
        D.ImageContainsSL2 E l := by
  obtain ⟨C₇, hC₇pos, Exc, hExc, h37⟩ := lemma_3_7 D KV hKV ε hε
  obtain ⟨B, hB⟩ := D.faltingsHeight_bddBelow
  have hKpos : (0:ℝ) < (23040:ℝ) ^ ε := Real.rpow_pos_of_pos (by norm_num) ε
  refine ⟨C₇ * (23040:ℝ) ^ ε + |B| + 1, by positivity, Exc, hExc, ?_⟩
  intro E l hl hExcE hcond
  set E' := D.torsionExt E with hE'def
  have hcls : D.cls E' = D.cls E := D.cls_torsionExt E
  have hss : D.SemiStable E' := D.semiStable_torsionExt E
  have hExcE' : D.cls E' ∉ Exc := by rw [hcls]; exact hExcE
  -- ★★★★**`5 ≤ l` はどちらの枝からも出る**（`Check/GenEll/VeluQuotOKNeedsL5.lean`、第 1434）
  have hl5 : 5 ≤ l := by
    rcases hcond with ⟨hnum, -⟩ | ⟨-, -, hcop⟩
    · -- ☆枝 (a): 括弧の中は `≥ 1`、`d ≥ 1` なので `l ≥ 2304000`
      set d : ℝ := (D.degOfDefinition E : ℝ) with hddef
      set F : ℝ := D.faltingsHeight (D.cls E) with hFdef
      set P : ℝ := d ^ ε with hPdef
      set K : ℝ := (23040:ℝ) ^ ε with hKdef
      have hd1 : (1:ℝ) ≤ d := by
        rw [hddef]; exact_mod_cast D.degOfDefinition_pos E
      have hP1 : (1:ℝ) ≤ P := by
        rw [hPdef]; exact Real.one_le_rpow hd1 hε.le
      have hK0 : (0:ℝ) < K := by rw [hKdef]; exact hKpos
      have hBF : -|B| ≤ F := le_trans (neg_abs_le B) (hB (D.cls E))
      have hBnn : (0:ℝ) ≤ |B| := abs_nonneg B
      have hC0 : (0:ℝ) ≤ C₇ * K := by positivity
      have hbr : (1:ℝ) ≤ F + (C₇ * K + |B| + 1) * P := by nlinarith
      have hlarge : (2304000:ℝ) ≤ (l:ℝ) := by nlinarith
      have h5 : (5:ℝ) ≤ (l:ℝ) := by linarith
      exact_mod_cast h5
    · -- ☆枝 (b): `30` と素な素数は `7` 以上
      by_contra hlt
      rw [not_le] at hlt
      interval_cases l
      · exact absurd hl (by decide)
      · exact absurd hl (by decide)
      · exact absurd hcop (by decide)
      · exact absurd hcop (by decide)
      · exact absurd hl (by decide)
  obtain ⟨ha, hb, hc⟩ := h37 E' l hl hl5 hss
    (100 * (D.degOfDefinition E' : ℝ)
        * (D.faltingsHeight (D.cls E') + C₇ * (D.degOfDefinition E' : ℝ) ^ ε) ≤ (l : ℝ)
      ∧ D.HasMultRed E')
    (D.cls E' ∈ KV ∧ D.PrimeToLocalHeights E' l) Iff.rfl Iff.rfl
  have hmain : D.HasMultRed E' ∧ D.PrimeToLocalHeights E' l ∧ ¬ D.HasLCyclic E' l := by
    rcases hcond with ⟨hnum, hpot⟩ | ⟨hKVE, hpl, hcop⟩
    · -- 条件 (a): 基底変換して `Lemma 3.7` の (a) を満たす
      have hmult' : D.HasMultRed E' := D.hasMultRed_torsionExt E hpot
      have hnumA : 100 * (D.degOfDefinition E' : ℝ)
          * (D.faltingsHeight (D.cls E') + C₇ * (D.degOfDefinition E' : ℝ) ^ ε) ≤ (l : ℝ) := by
        rw [hcls]
        set d : ℝ := (D.degOfDefinition E : ℝ) with hddef
        set d' : ℝ := (D.degOfDefinition E' : ℝ) with hd'def
        set F : ℝ := D.faltingsHeight (D.cls E) with hFdef
        have hd1 : (1:ℝ) ≤ d := by
          rw [hddef]; exact_mod_cast D.degOfDefinition_pos E
        have hd'1 : (1:ℝ) ≤ d' := by
          rw [hd'def]; exact_mod_cast D.degOfDefinition_pos E'
        have hdd' : d' ≤ 23040 * d := by
          rw [hd'def, hddef]
          exact_mod_cast D.degOfDefinition_torsionExt E
        have hdpos : (0:ℝ) < d := by linarith
        have hd'pos : (0:ℝ) < d' := by linarith
        set K : ℝ := (23040:ℝ) ^ ε with hKdef
        set P : ℝ := d ^ ε with hPdef
        set P' : ℝ := d' ^ ε with hP'def
        have hP1 : (1:ℝ) ≤ P := by rw [hPdef]; exact Real.one_le_rpow hd1 hε.le
        have hP'1 : (1:ℝ) ≤ P' := by rw [hP'def]; exact Real.one_le_rpow hd'1 hε.le
        have hK1 : (1:ℝ) ≤ K := by rw [hKdef]; exact Real.one_le_rpow (by norm_num) hε.le
        have hPKP : P' ≤ K * P := by
          rw [hP'def, hKdef, hPdef, ← Real.mul_rpow (by norm_num) hdpos.le]
          exact Real.rpow_le_rpow hd'pos.le hdd' hε.le
        have hBF : -|B| ≤ F := le_trans (neg_abs_le B) (hB (D.cls E))
        have hBnn : (0:ℝ) ≤ |B| := abs_nonneg B
        have hdPnn : (0:ℝ) ≤ d * P := by nlinarith
        have hBdPnn : (0:ℝ) ≤ |B| * (d * P) := mul_nonneg hBnn hdPnn
        have hstep : 100 * d' * (F + C₇ * P')
            ≤ 23040 * 100 * d * (F + (C₇ * K + |B| + 1) * P) := by
          -- ★両辺を単項へ展開しておく(`linarith` は積を原子として見る)
          have hLexp : 100 * d' * (F + C₇ * P') = 100 * d' * F + 100 * C₇ * (d' * P') := by ring
          have hRexp : 23040 * 100 * d * (F + (C₇ * K + |B| + 1) * P)
              = 23040 * 100 * d * F + 23040 * 100 * C₇ * K * (d * P)
                + 23040 * 100 * (|B| * (d * P)) + 23040 * 100 * (d * P) := by ring
          rw [hLexp, hRexp]
          -- ★第 2 項: `d'·d'^ε ≤ (23040d)·(23040^ε d^ε)`
          have hterm2 : 100 * C₇ * (d' * P') ≤ 23040 * 100 * C₇ * K * (d * P) := by
            have h1 : d' * P' ≤ (23040 * d) * (K * P) :=
              mul_le_mul hdd' hPKP (by linarith) (by linarith)
            have h2 := mul_le_mul_of_nonneg_left h1 (by positivity : (0:ℝ) ≤ 100 * C₇)
            linarith [h2]
          rcases le_or_gt 0 F with hF | hF
          · -- ★`F ≥ 0` なら第 1 項は `d' ≤ 23040d` で直ちに出る
            have hterm1 : 100 * d' * F ≤ 23040 * 100 * d * F := by nlinarith
            linarith [hterm1, hterm2, hBdPnn, hdPnn]
          · -- ★★`F < 0` なら第 1 項は逆向きになるが、超過分は `|B|·23040d` で押さえられ、
            -- `C` の `|B| + 1` の分がそれを吸収する(`P ≥ 1` が効く)
            have hFB : -F ≤ |B| := by linarith
            have hmul : (-F) * (23040 * d - d') ≤ |B| * (23040 * d) :=
              mul_le_mul hFB (by linarith) (by linarith) hBnn
            have hterm1 : 100 * d' * F - 23040 * 100 * d * F ≤ 23040 * 100 * (|B| * d) := by
              nlinarith [hmul]
            have hdP : d ≤ d * P := by nlinarith
            have hBdP : |B| * d ≤ |B| * (d * P) := mul_le_mul_of_nonneg_left hdP hBnn
            linarith [hterm1, hterm2, hBdP, hdPnn]
        linarith
      have hcA : (100 * (D.degOfDefinition E' : ℝ)
          * (D.faltingsHeight (D.cls E') + C₇ * (D.degOfDefinition E' : ℝ) ^ ε) ≤ (l : ℝ)
        ∧ D.HasMultRed E') := ⟨hnumA, hmult'⟩
      refine ⟨hmult', ha hcA, fun hcyc => hExcE' (hc (Or.inl hcA) hcyc)⟩
    · -- 条件 (b): `30` と互いに素なので局所高さと素であることが `L′` へ移る
      have hpl' : D.PrimeToLocalHeights E' l := D.primeToLocalHeights_torsionExt E l hpl hcop
      have hcB : (D.cls E' ∈ KV ∧ D.PrimeToLocalHeights E' l) := by
        refine ⟨?_, hpl'⟩
        rw [hcls]; exact hKVE
      exact ⟨hb hcB hExcE', hpl', fun hcyc => hExcE' (hc (Or.inr hcB) hcyc)⟩
  -- ★`Lemma 3.1, (iv)` が要求する `5 ≤ l` は上で作った（どちらの条件からも出る）
  exact D.imageContainsSL2_of_torsionExt E l hl hl5 hmain.1 hmain.2.1 hmain.2.2

/-! ## ★出典の紐付け(`.src`) -/

def theorem_3_8.src : Source :=
  { paper := "GenEll", pdfPage := 19, item := "Theorem 3.8",
    sectionId := "genell-thm-3-8" }

def theorem_3_8.needs : List ProofObligation :=
  [ .citation "[ABC3]" "lemma_3_7(Skeleton/GenEll/Section3.lean)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.lemma_3_7") 1,
    .citation "[ABC3]" "EllModuliData.imageContainsSL2_of_torsionExt(Interface の欄)"
      (.inProject "ABC3" "ABC3.Interface.GenEll.EllModuliData") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1445）——場所を規約どおりに直した**。" ++
       "☆証明は `sorry` 0 であったが本体が `Skeleton/` にあった。" ++
       "`Lemma 4.1` / `Lemma 4.2` / `Corollary 4.3` / `Corollary 4.4`（第 1444）と" ++
       "同じ形（`Skeleton` は写して `Found` へ委譲）に揃えた。" ++
       "★中身は 1 文字も変えていない。" ++
       "★★ただし `EllModuliData` の欄は `Interface` で受けたままであり、" ++
       "「絶対的に証明された」ことを意味しない——欄を埋めるのは `EllModuliWitness.lean` の仕事で、" ++
       "そこには `j(E′)` の整性と Tate 加群の Galois 像が残っている。") 1 ]

end ABC3.Found.GenEll
