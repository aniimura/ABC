import ABC3.Found.GaloisRep.PrimeOverX

/-!
# Galois (G5) 第 135 ブロック —— **★★★★★★素イデアルの生成元**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★`P = (x − c, y − y₀)` が出た

第 134 で `x − c ∈ P` を得た。★本ブロックは**残りの生成元**を出し、
`P` がちょうど 2 元で生成されることを示す。

| 段 | 内容 |
|---|---|
| `z ≡ ±s` | `z² = Ψ₂Sq(x) ≡ Ψ₂Sq(c) = s²` と `P` が素であること(第 133 の分解) |
| `y ≡ y₀` | `z = 2y + a₁x + a₃` を解く(`2` が単元) |
| **`P = (x − c, y − y₀)`** | 座標環の元が `x`,`y` の多項式であること |

## ★★★★最後の段の機構——**零点定理も基底も要らない**

`AdjoinRoot.mk` が全射なので、任意の `α ∈ F[W]` は `α = mk(p)`(`p ∈ F[X][Y]`)。
★`Q := (x − c, y − y₀)` で割ると `x ≡ c`・`y ≡ y₀` だから、
**環準同型の一致(`Polynomial.ringHom_ext`)**により

    α ≡ p(c, y₀) ∈ F   (mod Q)

★★したがって `α − p(c,y₀) ∈ Q ⊆ P`。`α ∈ P` なら `p(c,y₀) ∈ P` だが、
`F` の 0 でない元は単元だから `p(c,y₀) = 0`、すなわち `α ∈ Q`。

★★★**`Q` が極大であることを示す必要すら無い**——`P ⊆ Q` が直接出る。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `exists_genZ_sub_mem` | ★★★`∃ s, s² = Ψ₂Sq(c) ∧ z − s ∈ P` |
| `quotient_algebraMap_poly` | ★`F[X]` の元は `mod I` で定数になる |
| `exists_const_sub_mem` | ★★★★**`α − v ∈ I`(`v ∈ F`)** |
| `prime_eq_span` | ★★★★★**`P = (x − c, y − y₀)`** |
| `exists_genY_sub_mem` | ★★★`y − y₀ ∈ P`(`y₀ = (s − a₁c − a₃)/2`) |
| `exists_prime_gens` | ★★★★★★**まとめ——0 でない素イデアルの生成元** |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve WeierstrassCurve.Affine Polynomial

variable {F : Type} [Field F]

/-- ★★★**`z ≡ ±s` (mod `P`)**——`s² = Ψ₂Sq(c)`。

★第 133 の分解 `(z − s)(z + s) = (x − c)·g(x)` の右辺が `P` に入るので、
`P` が素であることから左辺のどちらかが `P` に入る。 -/
theorem exists_genZ_sub_mem [IsAlgClosed F] (W : WeierstrassCurve.Affine F)
    (P : Ideal W.CoordinateRing) [hPp : P.IsPrime] {c : F}
    (hc : genX W - algebraMap F W.CoordinateRing c ∈ P) :
    ∃ s : F, s ^ 2 = W.Ψ₂Sq.eval c ∧ genZ W - algebraMap F W.CoordinateRing s ∈ P := by
  obtain ⟨s₀, hs₀⟩ :=
    IsAlgClosed.exists_pow_nat_eq (k := F) (W.Ψ₂Sq.eval c) (n := 2) (by norm_num)
  obtain ⟨g, hg⟩ := genZ_sub_factor W hs₀
  have hmem : (genZ W - algebraMap F W.CoordinateRing s₀)
      * (genZ W + algebraMap F W.CoordinateRing s₀) ∈ P := by
    rw [hg]; exact Ideal.mul_mem_right _ _ hc
  rcases hPp.mem_or_mem hmem with h | h
  · exact ⟨s₀, hs₀, h⟩
  · refine ⟨-s₀, by rw [neg_pow]; simpa using hs₀, ?_⟩
    rwa [map_neg, sub_neg_eq_add]

/-- ★`x ≡ c` (mod `I`) なら `F[X]` の元は `mod I` で定数 `q(c)` になる。 -/
theorem quotient_algebraMap_poly (W : WeierstrassCurve.Affine F)
    (I : Ideal W.CoordinateRing) {c : F}
    (hx : genX W - algebraMap F W.CoordinateRing c ∈ I) (q : Polynomial F) :
    Ideal.Quotient.mk I (algebraMap (Polynomial F) W.CoordinateRing q)
      = Ideal.Quotient.mk I (algebraMap F W.CoordinateRing (q.eval c)) := by
  have hext : ((Ideal.Quotient.mk I).comp (algebraMap (Polynomial F) W.CoordinateRing))
      = (((Ideal.Quotient.mk I).comp (algebraMap F W.CoordinateRing)).comp
          (Polynomial.evalRingHom c)) := by
    refine Polynomial.ringHom_ext (fun a => ?_) ?_
    · show Ideal.Quotient.mk I (algebraMap (Polynomial F) W.CoordinateRing (C a))
        = Ideal.Quotient.mk I (algebraMap F W.CoordinateRing ((C a).eval c))
      rw [Polynomial.eval_C]
      rfl
    · show Ideal.Quotient.mk I (algebraMap (Polynomial F) W.CoordinateRing X)
        = Ideal.Quotient.mk I (algebraMap F W.CoordinateRing (X.eval c))
      rw [Polynomial.eval_X, Ideal.Quotient.eq]
      exact hx
  exact congrArg (fun f : Polynomial F →+* _ => f q) hext

/-- ★★★★**座標環の元は `mod (x − c, y − y₀)` で定数になる**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★`AdjoinRoot.mk` の全射性と `Polynomial.ringHom_ext` だけで出る。
★★**零点定理も、`F[X]` 上の基底も要らない。** -/
theorem exists_const_sub_mem (W : WeierstrassCurve.Affine F)
    (I : Ideal W.CoordinateRing) {c y₀ : F}
    (hx : genX W - algebraMap F W.CoordinateRing c ∈ I)
    (hy : genY W - algebraMap F W.CoordinateRing y₀ ∈ I) (α : W.CoordinateRing) :
    ∃ v : F, α - algebraMap F W.CoordinateRing v ∈ I := by
  obtain ⟨p, rfl⟩ : ∃ p, CoordinateRing.mk W p = α := AdjoinRoot.mk_surjective α
  refine ⟨Polynomial.eval₂ (Polynomial.evalRingHom c) y₀ p, ?_⟩
  rw [← Ideal.Quotient.eq]
  have hext : ((Ideal.Quotient.mk I).comp (CoordinateRing.mk W))
      = (((Ideal.Quotient.mk I).comp (algebraMap F W.CoordinateRing)).comp
          (Polynomial.eval₂RingHom (Polynomial.evalRingHom c) y₀)) := by
    refine Polynomial.ringHom_ext (fun q => ?_) ?_
    · show Ideal.Quotient.mk I (CoordinateRing.mk W (C q))
        = Ideal.Quotient.mk I (algebraMap F W.CoordinateRing
            (Polynomial.eval₂ (Polynomial.evalRingHom c) y₀ (C q)))
      rw [Polynomial.eval₂_C]
      exact quotient_algebraMap_poly W I hx q
    · show Ideal.Quotient.mk I (CoordinateRing.mk W X)
        = Ideal.Quotient.mk I (algebraMap F W.CoordinateRing
            (Polynomial.eval₂ (Polynomial.evalRingHom c) y₀ X))
      rw [Polynomial.eval₂_X, Ideal.Quotient.eq]
      exact hy
  exact congrArg (fun f : Polynomial (Polynomial F) →+* _ => f p) hext

/-- ★★★★★**素イデアルはちょうど `(x − c, y − y₀)` である**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★`Q` が極大であることを示す必要は無い——`α ∈ P` に対し `α ≡ v ∈ F` で
`v ∈ P` かつ `F` の 0 でない元は単元だから `v = 0`、すなわち `α ∈ Q`。 -/
theorem prime_eq_span (W : WeierstrassCurve.Affine F)
    (P : Ideal W.CoordinateRing) [hPp : P.IsPrime] {c y₀ : F}
    (hx : genX W - algebraMap F W.CoordinateRing c ∈ P)
    (hy : genY W - algebraMap F W.CoordinateRing y₀ ∈ P) :
    P = Ideal.span {genX W - algebraMap F W.CoordinateRing c,
      genY W - algebraMap F W.CoordinateRing y₀} := by
  set Q : Ideal W.CoordinateRing := Ideal.span {genX W - algebraMap F W.CoordinateRing c,
    genY W - algebraMap F W.CoordinateRing y₀} with hQ
  have hxQ : genX W - algebraMap F W.CoordinateRing c ∈ Q := Ideal.subset_span (by simp)
  have hyQ : genY W - algebraMap F W.CoordinateRing y₀ ∈ Q := Ideal.subset_span (by simp)
  have hQP : Q ≤ P := by
    rw [hQ, Ideal.span_le]
    rintro z (rfl | rfl) <;> simpa using ‹_›
  refine le_antisymm (fun α hα => ?_) hQP
  obtain ⟨v, hv⟩ := exists_const_sub_mem W Q hxQ hyQ α
  have hvP : algebraMap F W.CoordinateRing v ∈ P := by
    have : α - (α - algebraMap F W.CoordinateRing v) ∈ P := P.sub_mem hα (hQP hv)
    simpa using this
  have hv0 : v = 0 := by
    by_contra hne
    exact hPp.ne_top (Ideal.eq_top_of_isUnit_mem P hvP
      ((isUnit_iff_ne_zero.2 hne).map (algebraMap F W.CoordinateRing)))
  rw [hv0, map_zero, sub_zero] at hv
  exact hv

/-- ★★★**`y ≡ y₀ = (s − a₁c − a₃)/2` (mod `P`)**——`z = 2y + a₁x + a₃` を解く。 -/
theorem exists_genY_sub_mem [IsAlgClosed F] (h2 : IsUnit (2 : F))
    (W : WeierstrassCurve.Affine F) (P : Ideal W.CoordinateRing) [hPp : P.IsPrime] {c : F}
    (hc : genX W - algebraMap F W.CoordinateRing c ∈ P) :
    ∃ s : F, s ^ 2 = W.Ψ₂Sq.eval c ∧
      genZ W - algebraMap F W.CoordinateRing s ∈ P ∧
      genY W - algebraMap F W.CoordinateRing ((s - W.a₁ * c - W.a₃) / 2) ∈ P := by
  obtain ⟨s, hs, hz⟩ := exists_genZ_sub_mem W P hc
  refine ⟨s, hs, hz, ?_⟩
  have h2' : (2 : F) ≠ 0 := h2.ne_zero
  set A := algebraMap F W.CoordinateRing with hA
  set y₀ : F := (s - W.a₁ * c - W.a₃) / 2 with hy₀
  have hy0 : (2 : F) * y₀ = s - W.a₁ * c - W.a₃ := by
    rw [hy₀]; field_simp
  have hkey : (2 : W.CoordinateRing) * A y₀ = A s - A W.a₁ * A c - A W.a₃ := by
    rw [show (2 : W.CoordinateRing) = A 2 from (map_ofNat _ _).symm, ← map_mul, hy0,
      map_sub, map_sub, map_mul]
  have hid : (2 : W.CoordinateRing) * (genY W - A y₀)
      = (genZ W - A s) - A W.a₁ * (genX W - A c) := by
    rw [mul_sub, hkey, genZ]; ring
  have hmul : (2 : W.CoordinateRing) * (genY W - A y₀) ∈ P := by
    rw [hid]; exact P.sub_mem hz (Ideal.mul_mem_left _ _ hc)
  have hu : IsUnit (2 : W.CoordinateRing) := by
    rw [show (2 : W.CoordinateRing) = A 2 from (map_ofNat _ _).symm]
    exact h2.map A
  have := Ideal.mul_mem_left P (↑hu.unit⁻¹) hmul
  rwa [← mul_assoc, IsUnit.val_inv_mul, one_mul] at this

/-- ★★★★★★**0 でない素イデアルの生成元**——`P = (x − c, y − y₀)`。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★これが DVR 判定(第 132・133 の局所構造)の入力である。 -/
theorem exists_prime_gens [IsAlgClosed F] (h2 : IsUnit (2 : F))
    (W : WeierstrassCurve.Affine F) (P : Ideal W.CoordinateRing) (hP : P ≠ ⊥) [hPp : P.IsPrime] :
    ∃ c s : F, s ^ 2 = W.Ψ₂Sq.eval c ∧
      genZ W - algebraMap F W.CoordinateRing s ∈ P ∧
      P = Ideal.span {genX W - algebraMap F W.CoordinateRing c,
        genY W - algebraMap F W.CoordinateRing ((s - W.a₁ * c - W.a₃) / 2)} := by
  obtain ⟨c, hc⟩ := exists_genX_sub_mem W P hP
  obtain ⟨s, hs, hz, hy⟩ := exists_genY_sub_mem h2 W P hc
  exact ⟨c, s, hs, hz, prime_eq_span W P hc hy⟩

/-! ## ★出典の紐付け(`.src`) -/

def exists_prime_gens.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の構成——0 でない素イデアルが (x−c, y−y₀) であること)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
