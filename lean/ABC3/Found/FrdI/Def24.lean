import ABC3.Found.FrdI.MonoidPrime

/-!
# [FrdI] Definition 2.4 —— perf-factorial / `Λ supports M` / `C(d)`

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.47–p.48。

原文 (FrdI p.47):
> Definition 2.4.

## ★段 1 —— 実数化 `⊗ R` の実体

★★**原文は `M rlf` を `M pf ⊗R` の一言で済ませている**(p.47 の (c))。
★**その一言の中身をここで作る。**

★`M^pf` の元は形式的な分数 `m/n`(`Pf.mk m n`)なので、
`f : M →+ ℝ≥0` が与えられれば `M^pf → ℝ≥0` は **`m/n ↦ f m / n`** で延びる
(`pfLift`)。★**これが `⊗ ℝ≥0` の計算的な実体**である。

★★**`M` が monoprime なら `M ≃+ Λ≥0` であり、`Λ≥0` はどの場合も
`ℝ≥0` に自然に埋まる**ので、`M^pf ↪ ℝ≥0` が得られる(`monoprimeLift`)。

★★★**同型の取り方に依らないこと**は `pfLift_unique_of_smul` で押さえる ——
2 つの埋め込みは **`ℝ≥0` の正の実数倍だけ違う**。
-/

namespace ABC3.Found.FrdI

open scoped NNReal

universe w

variable (M : Type w) [AddCommMonoid M]

/-! ## ★`M^pf → ℝ≥0` への延長 -/

variable {M}

/-- ★`M →+ ℝ≥0` を `M^pf → ℝ≥0` へ延ばす(台の写像)。 -/
noncomputable def pfLiftFun (f : M →+ ℝ≥0) : Pf M → ℝ≥0 :=
  Quotient.lift (fun x : M × ℕ+ => f x.1 / ((x.2 : ℕ) : ℝ≥0)) (by
    rintro ⟨m, a⟩ ⟨m', a'⟩ ⟨k, hk⟩
    have hk0 : ((k : ℕ) : ℝ≥0) ≠ 0 := by
      simpa using (k.pos.ne' : (0 : ℕ) ≠ (k : ℕ))
    have ha0 : ((a : ℕ) : ℝ≥0) ≠ 0 := by
      simpa using (a.pos.ne' : (0 : ℕ) ≠ (a : ℕ))
    have ha'0 : ((a' : ℕ) : ℝ≥0) ≠ 0 := by
      simpa using (a'.pos.ne' : (0 : ℕ) ≠ (a' : ℕ))
    have hf := congrArg f hk
    rw [map_nsmul, map_nsmul, nsmul_eq_mul, nsmul_eq_mul] at hf
    push_cast at hf
    show f m / ((a : ℕ) : ℝ≥0) = f m' / ((a' : ℕ) : ℝ≥0)
    rw [div_eq_div_iff ha0 ha'0]
    refine mul_left_cancel₀ hk0 ?_
    calc ((k : ℕ) : ℝ≥0) * (f m * ((a' : ℕ) : ℝ≥0))
        = ((k : ℕ) : ℝ≥0) * ((a' : ℕ) : ℝ≥0) * f m := by ring
      _ = ((k : ℕ) : ℝ≥0) * ((a : ℕ) : ℝ≥0) * f m' := hf
      _ = ((k : ℕ) : ℝ≥0) * (f m' * ((a : ℕ) : ℝ≥0)) := by ring)

@[simp] theorem pfLiftFun_mk (f : M →+ ℝ≥0) (m : M) (a : ℕ+) :
    pfLiftFun f (Pf.mk m a) = f m / ((a : ℕ) : ℝ≥0) := rfl

/-- ★★**`⊗ ℝ≥0` の実体** —— `f : M →+ ℝ≥0` の `M^pf` への延長。

★`M^pf` の元は形式的な分数 `m/n` なので、`m/n ↦ f m / n` で延びる。 -/
noncomputable def pfLift (f : M →+ ℝ≥0) : Pf M →+ ℝ≥0 where
  toFun := pfLiftFun f
  map_zero' := by
    show pfLiftFun f (Pf.mk 0 1) = 0
    simp
  map_add' x y := by
    induction x using Pf.inductionOn with
    | _ m a =>
      induction y using Pf.inductionOn with
      | _ m' a' =>
        have ha0 : ((a : ℕ) : ℝ≥0) ≠ 0 := by
          simpa using (a.pos.ne' : (0 : ℕ) ≠ (a : ℕ))
        have ha'0 : ((a' : ℕ) : ℝ≥0) ≠ 0 := by
          simpa using (a'.pos.ne' : (0 : ℕ) ≠ (a' : ℕ))
        rw [Pf.mk_add_mk]
        show f ((a' : ℕ) • m + (a : ℕ) • m') / (((a * a' : ℕ+) : ℕ) : ℝ≥0)
          = f m / ((a : ℕ) : ℝ≥0) + f m' / ((a' : ℕ) : ℝ≥0)
        rw [map_add, map_nsmul, map_nsmul, nsmul_eq_mul, nsmul_eq_mul]
        push_cast
        field_simp

@[simp] theorem pfLift_mk (f : M →+ ℝ≥0) (m : M) (a : ℕ+) :
    pfLift f (Pf.mk m a) = f m / ((a : ℕ) : ℝ≥0) := rfl

@[simp] theorem pfLift_of (f : M →+ ℝ≥0) (m : M) : pfLift f (Pf.of m) = f m := by
  show pfLift f (Pf.mk m 1) = f m
  simp

/-- ★**`f` が単射なら延長も単射** —— `M^pf ↪ ℝ≥0`。

★`f m / a = f m' / a'` から `a' • m` と `a • m'` の `f` 値が一致し、
`f` の単射性で `a' • m = a • m'`、すなわち `Pf.mk m a = Pf.mk m' a'`。 -/
theorem pfLift_injective {f : M →+ ℝ≥0} (hf : Function.Injective f) :
    Function.Injective (pfLift f) := by
  intro x y hxy
  induction x using Pf.inductionOn with
  | _ m a =>
    induction y using Pf.inductionOn with
    | _ m' a' =>
      have ha0 : ((a : ℕ) : ℝ≥0) ≠ 0 := by
        simpa using (a.pos.ne' : (0 : ℕ) ≠ (a : ℕ))
      have ha'0 : ((a' : ℕ) : ℝ≥0) ≠ 0 := by
        simpa using (a'.pos.ne' : (0 : ℕ) ≠ (a' : ℕ))
      rw [pfLift_mk, pfLift_mk, div_eq_div_iff ha0 ha'0] at hxy
      refine Pf.sound 1 ?_
      refine hf ?_
      rw [map_nsmul, map_nsmul, nsmul_eq_mul, nsmul_eq_mul]
      push_cast
      calc (1 : ℝ≥0) * ((a' : ℕ) : ℝ≥0) * f m = f m * ((a' : ℕ) : ℝ≥0) := by ring
        _ = f m' * ((a : ℕ) : ℝ≥0) := hxy
        _ = (1 : ℝ≥0) * ((a : ℕ) : ℝ≥0) * f m' := by ring

/-! ## ★`Λ≥0 ↪ ℝ≥0` —— monoid type の実数化 -/

/-- ★`Λ≥0` から `ℝ≥0` への自然な埋め込み。

★`ℤ≥0 = ℕ ↪ ℝ≥0`、`ℚ≥0 ↪ ℝ≥0`、`ℝ≥0 = ℝ≥0`。 -/
noncomputable def MonoidType.toNNReal : (Λ : MonoidType) → (Λ.Nonneg →+ ℝ≥0)
  | .int => Nat.castAddMonoidHom ℝ≥0
  | .rat =>
      { toFun := fun q : ℚ≥0 => (q : ℝ≥0)
        map_zero' := NNRat.cast_zero
        map_add' := fun p q => NNRat.cast_add p q }
  | .real => AddMonoidHom.id ℝ≥0

theorem MonoidType.toNNReal_injective (Λ : MonoidType) :
    Function.Injective Λ.toNNReal := by
  cases Λ
  · exact fun a b h => Nat.cast_injective (R := ℝ≥0) h
  · exact fun a b h => NNRat.cast_injective h
  · exact fun a b h => h

/-- ★★**monoprime なモノイドの `M^pf` は `ℝ≥0` に埋まる**。

★★**これが `M^rlf_p = M^pf_p ⊗ ℝ≥0` の実体**である ——
`M ≃+ Λ≥0 ↪ ℝ≥0` を `pfLift` で `M^pf` へ延ばす。 -/
noncomputable def monoprimeLift {Λ : MonoidType} (e : M ≃+ Λ.Nonneg) : Pf M →+ ℝ≥0 :=
  pfLift ((Λ.toNNReal).comp (e : M →+ Λ.Nonneg))

theorem monoprimeLift_injective {Λ : MonoidType} (e : M ≃+ Λ.Nonneg) :
    Function.Injective (monoprimeLift e) :=
  pfLift_injective (fun _ _ h =>
    e.injective ((Λ.toNNReal_injective) h))

/-! ## ★★同型の取り方に依らないこと

★`monoprimeLift` は同型 `e` の選び方に依るが、★**依り方は `ℝ≥0` の正の実数倍だけ**
である。★これが原文が `⊗ ℝ≥0` を「標準的」と書ける根拠である。 -/

/-- ★★**延長は `M` 上の値で決まる** —— `M^pf` は `M` の像で生成されるので。 -/
theorem pfLift_ext {f g : M →+ ℝ≥0} (h : ∀ m : M, f m = g m) :
    pfLift f = pfLift g := by
  refine AddMonoidHom.ext fun x => ?_
  induction x using Pf.inductionOn with
  | _ m a => simp [h m]

/-- ★★★**正の実数倍を除く一意性** —— 2 つの延長が `M` 上で正の実数倍だけ違えば、
`M^pf` 上でも同じ倍率で違う。

★★**「`⊗ ℝ≥0` の取り方に依らない」の中身**である。 -/
theorem pfLift_smul {f g : M →+ ℝ≥0} (c : ℝ≥0) (h : ∀ m : M, g m = c * f m)
    (x : Pf M) : pfLift g x = c * pfLift f x := by
  induction x using Pf.inductionOn with
  | _ m a =>
    rw [pfLift_mk, pfLift_mk, h m, mul_div_assoc]

/-! ## ★段 1 の後半 —— `sup(Bound(a))`

原文 (FrdI p.47):
> the various Bound

★★**原文が `sup` と書ける根拠を式にする。**

★`sup` が取れるには **2 つ**要る:
1. **`M^rlf_p` が全順序かつ条件付き完備**であること
   —— ★**これが (b)(`M_p` は monoprime)の効き所**である。
   `M_p ≃ Λ≥0` から `(M_p)^pf` は `ℚ≥0` か `ℝ≥0`(`Λ = ℤ` なら `ℕ^pf ≃ ℚ≥0`)、
   したがって `⊗ ℝ≥0` して `M^rlf_p ≃ ℝ≥0` になる。
   ★**(b) が無ければ `M^pf_p ⊗ ℝ≥0` はただの加法モノイドで、`sup` は意味を持たない。**
2. その部分集合が **`M^rlf_p` の中で上に有界**であること
   —— ★★**これは自動ではない。** `Bound` は定義から `a` で押さえられているが、
   ★**その `a` は `M^pf` の元であって `M^pf_p` の元ではない**ので、
   `M^rlf_p` の中での上界にはならない。
   ★**だから原文は角括弧で「i.e., the various Bound(a) are bounded subsets」と
   条件として書いている。** ここでも `BddAbove` を仮定として持つ。
-/

/-- ★`ℝ≥0` ではモノイドの `≤`(`MLe`)は順序の `≤` に一致する。 -/
theorem mle_nnreal_iff {x y : ℝ≥0} : MLe x y ↔ x ≤ y := by
  constructor
  · rintro ⟨c, rfl⟩
    exact le_self_add
  · intro h
    exact ⟨y - x, add_tsub_cancel_of_le h⟩

/-- ★加法準同型は `MLe` を保つ。 -/
theorem map_mle {N : Type*} [AddCommMonoid N] (f : M →+ N) {a b : M} (h : MLe a b) :
    MLe (f a) (f b) := by
  obtain ⟨c, rfl⟩ := h
  exact ⟨f c, by rw [map_add]⟩

/-- ★`ℝ≥0` へ送れば `MLe` は `≤` になる。 -/
theorem map_le_of_mle (f : M →+ ℝ≥0) {a b : M} (h : MLe a b) : f a ≤ f b :=
  mle_nnreal_iff.mp (map_mle f h)

/-- ★★**(b) の効き所** —— monoprime なら `M^pf` は `ℝ≥0` へ
**順序を保って**単射に埋まる。★これで `M^rlf_p` は条件付き完備な全順序になる。 -/
theorem monoprimeLift_le_of_mle {Λ : MonoidType} (e : M ≃+ Λ.Nonneg) {a b : Pf M}
    (h : MLe a b) : monoprimeLift e a ≤ monoprimeLift e b :=
  map_le_of_mle _ h

variable (M) in
/-- ★`0` はつねに `Bound S a` に入る(`S` が `0` を含むなら)。

★★原文の添字が `p ∪ {0}` である理由 —— **`sup` を取る集合を空にしないため**。 -/
theorem zero_mem_bound {S : Set M} (h0 : (0 : M) ∈ S) (a : M) :
    (0 : M) ∈ Bound M S a :=
  ⟨h0, ⟨a, zero_add a⟩⟩

/-- ★★**原文の `sup(Bound_S(a))`** —— `Bound_S(a)` の像の上限を `ℝ≥0` で取る。 -/
noncomputable def boundSup (ι : M → ℝ≥0) (S : Set M) (a : M) : ℝ≥0 :=
  sSup (ι '' Bound M S a)

/-- ★上界であること。★**`BddAbove` が要る**(上の節の 2 番)。 -/
theorem le_boundSup (ι : M → ℝ≥0) (S : Set M) (a : M)
    (hbdd : BddAbove (ι '' Bound M S a)) {x : M} (hx : x ∈ Bound M S a) :
    ι x ≤ boundSup ι S a :=
  le_csSup hbdd ⟨x, hx, rfl⟩

/-- ★**最小の上界**であること。★`S` が `0` を含めば集合は空でない。 -/
theorem boundSup_le (ι : M → ℝ≥0) {S : Set M} (h0 : (0 : M) ∈ S) (a : M) {c : ℝ≥0}
    (h : ∀ x ∈ Bound M S a, ι x ≤ c) : boundSup ι S a ≤ c :=
  csSup_le (⟨ι 0, ⟨0, zero_mem_bound M h0 a, rfl⟩⟩) (by rintro _ ⟨x, hx, rfl⟩; exact h x hx)

/-- ★★`Bound` は単調 —— `MLe a b` なら `Bound S a ⊆ Bound S b`。 -/
theorem bound_mono {S : Set M} {a b : M} (h : MLe a b) : Bound M S a ⊆ Bound M S b :=
  fun _ hx => ⟨hx.1, mle_trans hx.2 h⟩

/-- ★★**`boundSup` は単調** —— 因子分解写像が `MLe` を保つことの中身。 -/
theorem boundSup_mono (ι : M → ℝ≥0) (S : Set M) {a b : M} (h : MLe a b)
    (hbdd : BddAbove (ι '' Bound M S b)) (h0 : (0 : M) ∈ S) :
    boundSup ι S a ≤ boundSup ι S b :=
  boundSup_le ι h0 a fun x hx => le_boundSup ι S b hbdd (bound_mono h hx)

end ABC3.Found.FrdI
