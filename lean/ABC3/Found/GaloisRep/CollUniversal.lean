import ABC3.Found.GaloisRep.TateCollinearTrunc

/-!
# Galois (G6) 第 251 ブロック —— **★★★★★★共線性の万有な環**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★★`ℤ[U,V,W]` を分母で局所化する

第 223 の `ℤ[A,W]` を 3 変数にしたもの。`q = U V W` で、3 点は
`(U, VW)`, `(V, UW)`, `(W, UV)`。

★分母になるのは **6 つの「底」**——`U, V, W, VW, UW, UV`——とその `Q^{m+1}` 倍である。

    collDenomBases : Fin 6 → CollBase := ![U, V, W, VW, UW, UV]
    collDenomSet := range (fun i => 1 − 底 i) ∪ range (fun (m,i) => 1 − Q^{m+1}·底 i)

★★`Fin 6` で添字づけたので、単元性の場合分けが **2 通り**で済む
(第 223 は 3 通りだったが、6 個を並べていたら 6 通りになっていた)。

## ★★特殊化に要る仮定

    ∀ i, IsUnit (1 − collPts u v w i)      (`collPts u v w = ![u, v, w, vw, uw, uv]`)

これは「3 点のどれも原点でなく、どの 2 点の和も原点でない」という条件である。
★Tate 一意化の準同型性を述べるときの自然な非退化条件にあたる。

## ★★★★★★★★帰着の形

    collDefect u v w q = 0
      ⟸ ∀ n, (kU kV kW)^n ∣ collDefectTrunc n kU kV kW (kU kV kW)

★残るのは**万有な環での整除性ひとつ**である。第 223 の
`tateDefect_eq_zero_of_univ` と同じ形になった。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `CollUniv` | ★★★★★共線性の万有な環 |
| `collEval`・`collSpecialize` | ★★★★★★特殊化 |
| `collUnits_k` | ★★万有な環の側では 12 個の単元条件が自動で成り立つ |
| `collDefect_eq_zero_of_univ` | ★★★★★★★★**万有な環での整除性ひとつに帰着** |
-/

namespace ABC3.Found.GaloisRep

open Complex Real MvPolynomial

/-! ## ★★★★★万有な環 -/

/-- ★共線性の万有な多項式環 `ℤ[U,V,W]`。 -/
abbrev CollBase : Type := MvPolynomial (Fin 3) ℤ

noncomputable def cU : CollBase := MvPolynomial.X 0
noncomputable def cV : CollBase := MvPolynomial.X 1
noncomputable def cW : CollBase := MvPolynomial.X 2

/-- ★`q = U V W`。 -/
noncomputable def cQ : CollBase := cU * cV * cW

/-- ★分母の「底」——3 点とその相方。 -/
noncomputable def collDenomBases : Fin 6 → CollBase := ![cU, cV, cW, cV * cW, cU * cW, cU * cV]

/-- ★分母となる元の集合。★`Fin 6` で添字づけたので場合分けが 2 通りで済む。 -/
noncomputable def collDenomSet : Set CollBase :=
  (Set.range fun i : Fin 6 => 1 - collDenomBases i)
    ∪ (Set.range fun p : ℕ × Fin 6 => 1 - cQ ^ (p.1 + 1) * collDenomBases p.2)

noncomputable def collDenoms : Submonoid CollBase := Submonoid.closure collDenomSet

/-- ★★★★★**共線性の万有な環**——`ℤ[U,V,W]` を分母で局所化したもの。 -/
abbrev CollUniv : Type := Localization collDenoms

/-! ## ★★★★★★特殊化 -/

variable {R : Type} [CommRing R] {I : Ideal R}

noncomputable def collEval (u v w : R) : CollBase →+* R :=
  MvPolynomial.eval₂Hom (Int.castRingHom R) ![u, v, w]

@[simp] theorem collEval_U (u v w : R) : collEval u v w cU = u := by
  rw [collEval, cU, MvPolynomial.eval₂Hom_X']; rfl

@[simp] theorem collEval_V (u v w : R) : collEval u v w cV = v := by
  rw [collEval, cV, MvPolynomial.eval₂Hom_X']; rfl

@[simp] theorem collEval_W (u v w : R) : collEval u v w cW = w := by
  rw [collEval, cW, MvPolynomial.eval₂Hom_X']; rfl

@[simp] theorem collEval_Q (u v w : R) : collEval u v w cQ = u * v * w := by
  rw [cQ, map_mul, map_mul, collEval_U, collEval_V, collEval_W]

/-- ★3 点とその相方(`R` の側)。 -/
noncomputable def collPts (u v w : R) : Fin 6 → R := ![u, v, w, v * w, u * w, u * v]

theorem collEval_denomBases (u v w : R) (i : Fin 6) :
    collEval u v w (collDenomBases i) = collPts u v w i := by
  fin_cases i <;> simp [collDenomBases, collPts]

theorem collEval_isUnit_denomSet [IsAdicComplete I R] (u v w q : R) (hq : q ∈ I)
    (huvw : u * v * w = q) (hcp : ∀ i, IsUnit (1 - collPts u v w i))
    {x : CollBase} (hx : x ∈ collDenomSet) : IsUnit (collEval u v w x) := by
  rcases hx with ⟨i, rfl⟩ | ⟨⟨m, i⟩, rfl⟩
  · rw [map_sub, map_one, collEval_denomBases]
    exact hcp i
  · rw [map_sub, map_one, map_mul, map_pow, collEval_Q, collEval_denomBases, huvw]
    exact isUnit_one_sub (I := I) (pow_succ_mul_mem hq m)

theorem collEval_isUnit_denoms [IsAdicComplete I R] (u v w q : R) (hq : q ∈ I)
    (huvw : u * v * w = q) (hcp : ∀ i, IsUnit (1 - collPts u v w i))
    (y : collDenoms) : IsUnit (collEval u v w (y : CollBase)) := by
  obtain ⟨x, hx⟩ := y
  refine Submonoid.closure_induction (s := collDenomSet)
    (motive := fun z _ => IsUnit (collEval u v w z)) ?_ ?_ ?_ hx
  · intro z hz
    exact collEval_isUnit_denomSet u v w q hq huvw hcp hz
  · simp
  · intro z₁ z₂ _ _ h₁ h₂
    rw [map_mul]
    exact h₁.mul h₂

/-- ★★★★★★**万有な環から任意の完備環への特殊化**。 -/
noncomputable def collSpecialize [IsAdicComplete I R] (u v w q : R) (hq : q ∈ I)
    (huvw : u * v * w = q) (hcp : ∀ i, IsUnit (1 - collPts u v w i)) : CollUniv →+* R :=
  IsLocalization.lift (collEval_isUnit_denoms (I := I) u v w q hq huvw hcp)

theorem collSpecialize_base [IsAdicComplete I R] (u v w q : R) (hq : q ∈ I)
    (huvw : u * v * w = q) (hcp : ∀ i, IsUnit (1 - collPts u v w i)) (x : CollBase) :
    collSpecialize u v w q hq huvw hcp (algebraMap CollBase CollUniv x) = collEval u v w x :=
  IsLocalization.lift_eq _ x

end ABC3.Found.GaloisRep

/-! ## ★万有な環の側の元と単元性 -/

namespace ABC3.Found.GaloisRep

open Complex Real MvPolynomial

noncomputable def kU : CollUniv := algebraMap CollBase CollUniv cU
noncomputable def kV : CollUniv := algebraMap CollBase CollUniv cV
noncomputable def kW : CollUniv := algebraMap CollBase CollUniv cW

theorem kQ_eq : kU * kV * kW = algebraMap CollBase CollUniv cQ := by
  rw [cQ, map_mul, map_mul, kU, kV, kW]

theorem mem_collDenoms_of_mem (x : CollBase) (hx : x ∈ collDenomSet) : x ∈ collDenoms :=
  Submonoid.subset_closure hx

theorem isUnit_one_sub_kbase (i : Fin 6) :
    IsUnit (1 - algebraMap CollBase CollUniv (collDenomBases i)) := by
  have h : (1 : CollBase) - collDenomBases i ∈ collDenoms :=
    mem_collDenoms_of_mem _ (Or.inl ⟨i, rfl⟩)
  have h2 := IsLocalization.map_units (M := collDenoms) CollUniv (⟨_, h⟩ : collDenoms)
  simpa using h2

theorem isUnit_one_sub_kq (m : ℕ) (i : Fin 6) :
    IsUnit (1 - (algebraMap CollBase CollUniv cQ) ^ (m + 1)
      * algebraMap CollBase CollUniv (collDenomBases i)) := by
  have h : (1 : CollBase) - cQ ^ (m + 1) * collDenomBases i ∈ collDenoms :=
    mem_collDenoms_of_mem _ (Or.inr ⟨(m, i), rfl⟩)
  have h2 := IsLocalization.map_units (M := collDenoms) CollUniv (⟨_, h⟩ : collDenoms)
  simpa using h2

/-- ★★万有な環の側では 12 個の単元条件が自動で成り立つ。 -/
theorem collUnits_k : CollUnits kU kV kW := by
  have e0 : algebraMap CollBase CollUniv (collDenomBases 0) = kU := rfl
  have e1 : algebraMap CollBase CollUniv (collDenomBases 1) = kV := rfl
  have e2 : algebraMap CollBase CollUniv (collDenomBases 2) = kW := rfl
  have e3 : algebraMap CollBase CollUniv (collDenomBases 3) = kV * kW := by
    show algebraMap CollBase CollUniv (cV * cW) = kV * kW
    rw [map_mul]; rfl
  have e4 : algebraMap CollBase CollUniv (collDenomBases 4) = kU * kW := by
    show algebraMap CollBase CollUniv (cU * cW) = kU * kW
    rw [map_mul]; rfl
  have e5 : algebraMap CollBase CollUniv (collDenomBases 5) = kU * kV := by
    show algebraMap CollBase CollUniv (cU * cV) = kU * kV
    rw [map_mul]; rfl
  exact
    { hu := e0 ▸ isUnit_one_sub_kbase 0
      hv := e1 ▸ isUnit_one_sub_kbase 1
      hw := e2 ▸ isUnit_one_sub_kbase 2
      hvw := e3 ▸ isUnit_one_sub_kbase 3
      huw := e4 ▸ isUnit_one_sub_kbase 4
      huv := e5 ▸ isUnit_one_sub_kbase 5
      hqu := fun m => by rw [kQ_eq, ← e0]; exact isUnit_one_sub_kq m 0
      hqv := fun m => by rw [kQ_eq, ← e1]; exact isUnit_one_sub_kq m 1
      hqw := fun m => by rw [kQ_eq, ← e2]; exact isUnit_one_sub_kq m 2
      hqvw := fun m => by rw [kQ_eq, ← e3]; exact isUnit_one_sub_kq m 3
      hquw := fun m => by rw [kQ_eq, ← e4]; exact isUnit_one_sub_kq m 4
      hquv := fun m => by rw [kQ_eq, ← e5]; exact isUnit_one_sub_kq m 5 }

/-! ## ★★★★★★★★共線性の帰着 -/

variable {R : Type} [CommRing R] {I : Ideal R}

theorem collSpecialize_kU [IsAdicComplete I R] (u v w q : R) (hq : q ∈ I)
    (huvw : u * v * w = q) (hcp : ∀ i, IsUnit (1 - collPts u v w i)) :
    collSpecialize u v w q hq huvw hcp kU = u := by
  rw [kU, collSpecialize_base, collEval_U]

theorem collSpecialize_kV [IsAdicComplete I R] (u v w q : R) (hq : q ∈ I)
    (huvw : u * v * w = q) (hcp : ∀ i, IsUnit (1 - collPts u v w i)) :
    collSpecialize u v w q hq huvw hcp kV = v := by
  rw [kV, collSpecialize_base, collEval_V]

theorem collSpecialize_kW [IsAdicComplete I R] (u v w q : R) (hq : q ∈ I)
    (huvw : u * v * w = q) (hcp : ∀ i, IsUnit (1 - collPts u v w i)) :
    collSpecialize u v w q hq huvw hcp kW = w := by
  rw [kW, collSpecialize_base, collEval_W]

/-- ★★★★★★★★**共線性は万有な環での整除性ひとつに帰着した**。

`ℤ[U,V,W]` を 6 つの底とその `Q^{m+1}` 倍で局所化した環の中で

    (UVW)^n ∣ collDefectTrunc n U V W (UVW)

を示せば、**任意の完備環で 3 点は共線になる**。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem collDefect_eq_zero_of_univ [IsAdicComplete I R] (u v w q : R) (hq : q ∈ I)
    (huvw : u * v * w = q) (hcp : ∀ i, IsUnit (1 - collPts u v w i))
    (H : ∀ n : ℕ, ((kU * kV * kW) ^ n) ∣ collDefectTrunc n kU kV kW (kU * kV * kW)) :
    collDefect u v w q hq = 0 :=
  collDefect_eq_zero_of_specialize u v w q hq huvw kU kV kW
    (collSpecialize u v w q hq huvw hcp)
    (collSpecialize_kU u v w q hq huvw hcp) (collSpecialize_kV u v w q hq huvw hcp)
    (collSpecialize_kW u v w q hq huvw hcp) collUnits_k H

/-! ## ★出典の紐付け(`.src`) -/

def CollUniv.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——共線性の万有な環)",
    sectionId := "genell-def-3-3" }

def collDefect_eq_zero_of_univ.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——共線性を万有な環の整除性に帰着)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
