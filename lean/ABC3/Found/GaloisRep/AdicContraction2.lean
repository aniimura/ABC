import ABC3.Found.GaloisRep.TateInjUnit

/-!
# Galois (G6) 第 268 ブロック —— **★★★★★★2 変数の縮小写像定理**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★★★環帯の領域には 2 変数が要る

葉 (e) の**単元の領域**(第 266)では、`a` が単元なので相方 `w = q·a⁻¹` を
**割り算で作れた**——未知数は 1 つで済んだ。

★★環帯の領域(`a, w ∈ I`)では `a` が単元でないので `q/a` が作れない。
そこで `(a, w)` を**独立な 2 変数**として

    a = (x + y) − δ_S(a,w)      (`X + Y ≡ a`)
    w = −y + δ_Y(a,w)           (`Y ≡ −w`)

の不動点を取る。★★制約 `a·w = q` は**あとから**「差が `q` を決める」ことで回復する。

## ★★★2 変数版は 1 変数版の写し

第 102 の証明(反復 → Cauchy → `IsPrecomplete` → `IsHausdorff`)を対で走らせるだけ。
★成分ごとに `IsPrecomplete.prec'` を当てて極限 `A`, `B` を取り、
`F (f n).1 (f n).2 − F A B ∈ I^{n+1}` から不動点の式を出す。

★`nth_rewrite 2` で「引く側の `f (m+1)`」だけを書き換えるのが要点
(第 102 と同じ手)。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `exists_fixedPoint_of_contraction₂` | ★★★★★★**2 変数の縮小写像定理** |
-/

namespace ABC3.Found.GaloisRep

open Complex Real

variable {R : Type} [CommRing R] {I : Ideal R}

set_option maxHeartbeats 1200000 in
/-- ★★★★★★**2 変数の縮小写像定理**——第 102 を対で走らせたもの。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem exists_fixedPoint_of_contraction₂ [IsAdicComplete I R] (F G : R → R → R)
    (hFI : ∀ x ∈ I, ∀ y ∈ I, F x y ∈ I) (hGI : ∀ x ∈ I, ∀ y ∈ I, G x y ∈ I)
    (hconF : ∀ x ∈ I, ∀ y ∈ I, ∀ x' ∈ I, ∀ y' ∈ I, ∀ k : ℕ,
      x - x' ∈ I ^ k → y - y' ∈ I ^ k → F x y - F x' y' ∈ I ^ (k + 1))
    (hconG : ∀ x ∈ I, ∀ y ∈ I, ∀ x' ∈ I, ∀ y' ∈ I, ∀ k : ℕ,
      x - x' ∈ I ^ k → y - y' ∈ I ^ k → G x y - G x' y' ∈ I ^ (k + 1)) :
    ∃ a ∈ I, ∃ b ∈ I, F a b = a ∧ G a b = b := by
  set step : R × R → R × R := fun p => (F p.1 p.2, G p.1 p.2) with hstepdef
  set f : ℕ → R × R := fun n => step^[n] (0, 0) with hf
  have hsucc : ∀ n, f (n + 1) = step (f n) := by
    intro n
    rw [hf]
    simp only
    rw [Function.iterate_succ_apply' (f := step) (n := n)]
  have hA : ∀ n, (f (n + 1)).1 = F (f n).1 (f n).2 := fun n => by rw [hsucc n]
  have hB : ∀ n, (f (n + 1)).2 = G (f n).1 (f n).2 := fun n => by rw [hsucc n]
  have hfI : ∀ n, (f n).1 ∈ I ∧ (f n).2 ∈ I := by
    intro n
    induction n with
    | zero => rw [hf]; simp
    | succ m ih =>
      rw [hA m, hB m]
      exact ⟨hFI _ ih.1 _ ih.2, hGI _ ih.1 _ ih.2⟩
  have hstep : ∀ n : ℕ, (f (n + 1)).1 - (f n).1 ∈ I ^ n ∧ (f (n + 1)).2 - (f n).2 ∈ I ^ n := by
    intro n
    induction n with
    | zero => simp
    | succ m ih =>
      constructor
      · rw [hA (m + 1)]
        nth_rewrite 2 [hA m]
        exact hconF _ (hfI (m + 1)).1 _ (hfI (m + 1)).2 _ (hfI m).1 _ (hfI m).2 m ih.1 ih.2
      · rw [hB (m + 1)]
        nth_rewrite 2 [hB m]
        exact hconG _ (hfI (m + 1)).1 _ (hfI (m + 1)).2 _ (hfI m).1 _ (hfI m).2 m ih.1 ih.2
  have hcauchyA : ∀ {m n : ℕ}, m ≤ n →
      (f m).1 ≡ (f n).1 [SMOD (I ^ m • ⊤ : Submodule R R)] := by
    intro m n hmn
    induction n, hmn using Nat.le_induction with
    | base => rfl
    | succ n hmn ih =>
      refine ih.trans ?_
      rw [SModEq.sub_mem]
      have hle : I ^ n ≤ I ^ m := Ideal.pow_le_pow_right hmn
      simpa using neg_mem (hle (hstep n).1)
  have hcauchyB : ∀ {m n : ℕ}, m ≤ n →
      (f m).2 ≡ (f n).2 [SMOD (I ^ m • ⊤ : Submodule R R)] := by
    intro m n hmn
    induction n, hmn using Nat.le_induction with
    | base => rfl
    | succ n hmn ih =>
      refine ih.trans ?_
      rw [SModEq.sub_mem]
      have hle : I ^ n ≤ I ^ m := Ideal.pow_le_pow_right hmn
      simpa using neg_mem (hle (hstep n).2)
  obtain ⟨A, hAlim⟩ := IsPrecomplete.prec' (fun n => (f n).1) hcauchyA
  obtain ⟨B, hBlim⟩ := IsPrecomplete.prec' (fun n => (f n).2) hcauchyB
  have hAn : ∀ n, (f n).1 - A ∈ I ^ n := by
    intro n
    have hx := hAlim n
    rw [SModEq.sub_mem] at hx
    simpa using hx
  have hBn : ∀ n, (f n).2 - B ∈ I ^ n := by
    intro n
    have hx := hBlim n
    rw [SModEq.sub_mem] at hx
    simpa using hx
  have hAI : A ∈ I := by
    have h4 : A = (f 1).1 - ((f 1).1 - A) := by ring
    rw [h4]
    exact Submodule.sub_mem _ (hfI 1).1 (by simpa using hAn 1)
  have hBI : B ∈ I := by
    have h4 : B = (f 1).2 - ((f 1).2 - B) := by ring
    rw [h4]
    exact Submodule.sub_mem _ (hfI 1).2 (by simpa using hBn 1)
  refine ⟨A, hAI, B, hBI, ?_, ?_⟩
  · refine ((IsHausdorff.eq_iff_smodEq (I := I)).2 (fun n => ?_)).symm
    rw [SModEq.sub_mem]
    have hFdiff : F (f n).1 (f n).2 - F A B ∈ I ^ (n + 1) :=
      hconF _ (hfI n).1 _ (hfI n).2 _ hAI _ hBI n (hAn n) (hBn n)
    have hexp : A - F A B = -((f (n + 1)).1 - A) + (F (f n).1 (f n).2 - F A B) := by
      rw [hA n]; ring
    have hmem : A - F A B ∈ I ^ (n + 1) := by
      rw [hexp]
      exact Submodule.add_mem _ (neg_mem (hAn (n + 1))) hFdiff
    simpa using Ideal.pow_le_pow_right (Nat.le_succ n) hmem
  · refine ((IsHausdorff.eq_iff_smodEq (I := I)).2 (fun n => ?_)).symm
    rw [SModEq.sub_mem]
    have hGdiff : G (f n).1 (f n).2 - G A B ∈ I ^ (n + 1) :=
      hconG _ (hfI n).1 _ (hfI n).2 _ hAI _ hBI n (hAn n) (hBn n)
    have hexp : B - G A B = -((f (n + 1)).2 - B) + (G (f n).1 (f n).2 - G A B) := by
      rw [hB n]; ring
    have hmem : B - G A B ∈ I ^ (n + 1) := by
      rw [hexp]
      exact Submodule.add_mem _ (neg_mem (hBn (n + 1))) hGdiff
    simpa using Ideal.pow_le_pow_right (Nat.le_succ n) hmem

/-! ## ★出典の紐付け(`.src`) -/

def exists_fixedPoint_of_contraction₂.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——2 変数の縮小写像定理)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
