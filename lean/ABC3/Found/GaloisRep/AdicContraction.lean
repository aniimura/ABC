import ABC3.Found.GaloisRep.TateInvert

/-!
# Galois (G6) 第 102 ブロック —— **★★★★★完備 adic 環の縮小写像定理**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★★第 100 ブロックの道具を一般化する

Tate 母数の構成(第 100)では `q = t·v(q)` の不動点を作った。
★その機構は **`I` を `I^{k+1}` へ縮める写像**なら何でも動く:

    F : R → R,  x, y ∈ I,  x − y ∈ I^k  ⟹  F(x) − F(y) ∈ I^{k+1}
      ⟹  F は `I` の中に不動点を持つ

★★証明は `x_{n+1} = F(x_n)`(`x_0 = 0`)の反復であり、
連続する 2 項の差が `I^n` に入るので `IsPrecomplete` で極限が取れる。

★★★**位相を使わない**——`IsAdicComplete` だけで完結する。

## ★★これが次の層で繰り返し使われる

★Weierstrass 曲線の**形式群**を作るときの `w(z)`、
★★Hensel 型の引き上げ、★★★Tate 一意化の全射性(葉 (e))。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `exists_fixedPoint_of_contraction` | ★★★★★**縮小写像は不動点を持つ** |
-/

namespace ABC3.Found.GaloisRep

variable {R : Type} [CommRing R] {I : Ideal R}

/-- ★★★★★**完備 adic 環の縮小写像定理**。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

★第 100 ブロックの Tate 母数の構成を一般化したものである。 -/
theorem exists_fixedPoint_of_contraction [IsAdicComplete I R] (F : R → R)
    (hFI : ∀ x ∈ I, F x ∈ I) (h0 : F 0 ∈ I)
    (hcon : ∀ x ∈ I, ∀ y ∈ I, ∀ k : ℕ, x - y ∈ I ^ k → F x - F y ∈ I ^ (k + 1)) :
    ∃ q ∈ I, F q = q := by
  set f : ℕ → R := fun n => F^[n] 0 with hf
  have hsucc : ∀ n, f (n + 1) = F (f n) := by
    intro n
    rw [hf]
    simp only
    rw [Function.iterate_succ_apply' (f := F) (n := n)]
  have hfI : ∀ n, f n ∈ I := by
    intro n
    induction n with
    | zero =>
      rw [hf]
      simp
    | succ m ih =>
      rw [hsucc]
      exact hFI _ ih
  have hstep : ∀ n : ℕ, f (n + 1) - f n ∈ I ^ n := by
    intro n
    induction n with
    | zero => simp
    | succ m ih =>
      rw [hsucc (m + 1)]
      nth_rewrite 2 [hsucc m]
      exact hcon _ (hfI (m + 1)) _ (hfI m) m ih
  have hcauchy : ∀ {m n : ℕ}, m ≤ n → f m ≡ f n [SMOD (I ^ m • ⊤ : Submodule R R)] := by
    intro m n hmn
    induction n, hmn using Nat.le_induction with
    | base => rfl
    | succ n hmn ih =>
      refine ih.trans ?_
      rw [SModEq.sub_mem]
      have h := hstep n
      have hle : I ^ n ≤ I ^ m := Ideal.pow_le_pow_right hmn
      simpa using neg_mem (hle h)
  obtain ⟨L, hL⟩ := IsPrecomplete.prec' f hcauchy
  have hLI : L ∈ I := by
    have h1 := hL 1
    rw [SModEq.sub_mem] at h1
    have h2 : f 1 - L ∈ I := by simpa using h1
    have h4 : L = f 1 - (f 1 - L) := by ring
    rw [h4]
    exact Submodule.sub_mem _ (hfI 1) h2
  refine ⟨L, hLI, ?_⟩
  refine ((IsHausdorff.eq_iff_smodEq (I := I)).2 (fun n => ?_)).symm
  rw [SModEq.sub_mem]
  have hfn : f n - L ∈ I ^ n := by
    have hx := hL n
    rw [SModEq.sub_mem] at hx
    simpa using hx
  have hFdiff : F (f n) - F L ∈ I ^ (n + 1) := hcon _ (hfI n) _ hLI n hfn
  have hfn1 : f (n + 1) - L ∈ I ^ (n + 1) := by
    have hx := hL (n + 1)
    rw [SModEq.sub_mem] at hx
    simpa using hx
  have hexp : L - F L = -(f (n + 1) - L) + (F (f n) - F L) := by
    rw [hsucc]
    ring
  have hmem : L - F L ∈ I ^ (n + 1) := by
    rw [hexp]
    exact Submodule.add_mem _ (neg_mem hfn1) hFdiff
  have hle : I ^ (n + 1) ≤ I ^ n := Ideal.pow_le_pow_right (Nat.le_succ n)
  simpa using hle hmem

/-! ## ★出典の紐付け(`.src`) -/

def exists_fixedPoint_of_contraction.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(完備 adic 環の縮小写像定理)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
