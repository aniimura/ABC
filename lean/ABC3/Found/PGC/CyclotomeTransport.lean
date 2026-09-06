import ABC3.Found.PGC.TorsionCyclotome
import ABC3.Found.PGC.CyclotomicRecovery

/-!
# 円分子の移送 —— Proposition 1.1 を LCFT の 1 補題に還元する

[pGC] Proposition 1.1 への**経路 Λ**の節点 Λ10。Λ1
(`Found/PGC/TopAbelianization.lean`)と Λ2(`Found/PGC/TorsionCyclotome.lean`)の上に乗る。

到達点は

  `cyclotomicRecovery_of_torsionCyclotome :
     TorsionCyclotomeIsCyclotomic p → (cyclotomicCharacterObject (p := p)).RecoverableFromAbsGal`

である。すなわち **Proposition 1.1 は `TorsionCyclotomeIsCyclotomic p` という
名前の付いた 1 本の補題に還元された**。この 1 本は局所類体論(Lubin–Tate 相互律)が
与えるはずのもので、本ファイルでは証明していない。

## 筋

`α : Γ_K ≃ₜ* Γ_{K'}` と `g : Γ_K`、`n : ℕ` を取る。`m := p^n` と置く。

1. **両側で `μ_m` を止める開正規部分群を同時に取る**(Λ10 の本体)。
   `muFixer K m := {g : Γ_K | ∀ x, x^m = 1 → g x = x}`(`μ_m` を各点固定する元)は
   開正規部分群であり、`H := muFixer K m ⊓ α⁻¹(muFixer K' m)` も開正規。
   `H ≤ muFixer K m` と `α(H) ≤ muFixer K' m` が同時に成り立つ。
   ★これは `Found/PGC/DegreeTransport.lean` の `finrank_eq_of_absGal_equiv` が
   使ったのと**まったく同じ手**である——canonical でない開部分群で済ませる。
   違いは、あちらが `K(ζ)` の固定部分群(正規とは限らない)で足りたのに対し、
   こちらは共役作用を使うので**正規性が要る**点で、そのために
   `fixingSubgroup` ではなく `muFixer`(定義から正規)を使っている。
2. **円分子を移送する**。Λ2 の `cyclotomeEquiv (subgroupMapCME α H)` が
   `Λ_n(H) ≃* Λ_n(α H)` を与え、`cyclotomeEquiv_subgroupMapCME_conj` が
   `g` の作用を `α g` の作用に移すことを言う。
3. **仮定を両側に当てて指標を読む**。`TorsionCyclotomeIsCyclotomic` が
   `Λ_n(H) ≃* ℤ/p^n`(作用は `χ_{K,n}` 倍)と `Λ_n(α H) ≃* ℤ/p^n`(作用は `χ_{K',n}` 倍)
   を与えるので、合成 `Φ : ℤ/p^n ≃ ℤ/p^n` が
   `Φ(y^{χ_n(g)}) = (Φ y)^{χ'_n(α g)}` を満たす。`Φ` は群同型なので左辺は `(Φ y)^{χ_n(g)}`、
   よって `z^{χ_n(g)} = z^{χ'_n(α g)}` が全ての `z` について成り立ち、
   `z = ofAdd 1` を入れて `χ_n(g) = χ'_n(α g)`。
4. `n` を動かして `PadicInt.ext_of_toZModPow` で `ℤ_p` の等式に上げ、
   `cyclotomicCharacterObject_recoverable_iff` で結論。

## ★罠を守ったこと(2 つ)

* **`H` は結論の statement に現れない。** 結論は
  `(cyclotomicCharacterObject (p := p)).RecoverableFromAbsGal` であり、部分群は
  仮定(`TorsionCyclotomeIsCyclotomic` の `∀ S`)と証明の内側にしか出てこない。
  Proposition 1.2 の `∀ RD` で起きた「自由なパラメータが結論に現れる」退化
  (`Check/PGC/Prop12ForallRD.lean`、10 例目)を踏んでいない。
* **`Abelianization` を使っていない。** Λ1 が終始
  `TopologicalAbelianization`(交換子群の**位相的**閉包による商)を使う。

## ★核だけでは足りない、という別の確認との関係

`Found/PGC/CyclotomicRecovery.lean` の docstring は
「`ker(χ_m)` を回復しても `χ_m` は像の自己同型のぶんだけ不定」と記録していた。
本ファイルはその不定性を**円分子の生成元**で切っている:`Λ_n(H)` は
`Γ_K`-加群として `μ_{p^n}` と同型であり、同型は生成元の取り方に依存するが、
**`Γ` の作用の仕方**(すなわち指標)は同型の取り方に依らない。
これが `pow_eq_pow_of_equivariant` が形式化している内容である。
したがって本ファイルは「`ker(χ_n)` の群論的特徴づけ」も「`p` 冪の計数上界」も使わない。

## ★退化の自己検査——仮定を落とすと何が起きるか

* **`IsOpen S` を落とすと `TorsionCyclotomeIsCyclotomic` は偽になる**(自明化ではない)。
  `S = ⊥` は正規で `S ≤ muFixer F m` を(空虚に)満たすが、`↥⊥` は自明群なので
  `Λ_n(⊥)` も自明であり、`p^n > 1` のとき `ℤ/p^n` と同型になりえない
  (`not_nonempty_cyclotomeEquiv_bot`)。開性はこの退化をちょうど排除する。
* **`S ≤ muFixer F (p^n)` を落とすと偽になる**。`p` を奇素数、`F = ℚ_p`、`n = 1`、`S = ⊤` を取る。
  `μ_p ⊄ ℚ_p` なので `S ≤ muFixer F p` は成り立たない。局所類体論で
  `Γ_{ℚ_p}^{ab} ≅ (ℚ_p^×)^∧ ≅ ℤ^∧ × μ_{p−1} × ℤ_p` であり `p` 捩れは自明、
  したがって `Λ_1(⊤)` は自明で `ℤ/p` と同型にならない。
  (この計算は本ファイルでは形式化していない——`Found/PGC/UnitsPowP.lean` 系の在庫が
  `[𝒪^×:(𝒪^×)^p] = p^d · #μ_p` を持つので届く見込みだが、今は必要でない。)
* **`S.Normal` を落とすと型が付かない**。共役作用 `cyclotomeConj` が
  `H` の正規性を要求するからである(`Found/PGC/TorsionCyclotome.lean`)。
* **`n` を動かさないと結論に届かない**。`PadicInt.ext_of_toZModPow` は
  すべての `n` を要求する。ひとつの `n` だけでは `χ mod p^n` しか決まらない。

## 逸脱の記録

1. **`μ_{p^n}` の代わりに `Multiplicative (ℤ/p^n)` を使った。**
   規模測定の節点記述は「`tors_{p^n}(Γ_L^{ab}) ≅ μ_{p^n}(Aut(L/K)-同変)」だったが、
   `μ_{p^n} ⊆ L` のとき `μ_{p^n}(L)` は位数 `p^n` の巡回群であり、生成元を選べば
   `Multiplicative (ℤ/p^n)` と同型で、Galois 作用は `z ↦ z^{χ_n(g)}` になる。
   本ファイルはその形で書いた——`rootsOfUnity` の部分群としての配管
   (`Units.map` と `rootsOfUnity` の保存)を一段省ける。内容は同値である。
2. **`Aut(L/K)`-同変性は `Γ_K` の共役作用として書いた。**
   `Found/PGC/TorsionCyclotome.lean::cyclotomeConj_coe_self` が
   「`H` の元は自明に作用する」を与えるので、作用は `Γ_K ⧸ H = Gal(L_H/K)` を経由する。
   すなわち「`Aut(L/K)`-同変」と同じ内容である。底体 `K` を保持したのは、
   `L` だけを見ると `Aut(L/K)` を述べられないからである。
3. 原典(pGC §1)の論拠(Serre の局所類体論)を経由しない点は変わらない。
   `ResearchPaper/pgc-goal.md` に記録済みの逸脱の範囲内。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC ABC3.Found

variable {p : ℕ} [Fact p.Prime]

/-! ## 1. `μ_m` を各点固定する部分群 -/

/-- **`μ_m` を各点固定する `Γ_K` の部分群**。

`L_H ⊇ μ_m` を「`H` の元が `m` 乗根をすべて止める」と言い換えたもの
(`le_muFixer_iff`)。★`ker(χ mod m)` と同じものだが、**円分指標の連続性を
経由せずに**定義できるのが利点である(開性は `K(ζ)` の固定部分群を含むことから直接出る)。 -/
def muFixer (K : PAdicLocalField p) (m : ℕ) : Subgroup K.absGal where
  carrier := {g : K.absGal | ∀ x : K.closure, x ^ m = 1 → g x = x}
  one_mem' := fun _ _ => rfl
  mul_mem' := by
    intro a b ha hb x hx
    show a (b x) = x
    rw [hb x hx, ha x hx]
  inv_mem' := by
    intro a ha x hx
    have h1 : a x = x := ha x hx
    show a⁻¹ x = x
    conv_lhs => rw [← h1]
    show (a⁻¹ * a) x = x
    rw [inv_mul_cancel]; rfl

theorem mem_muFixer {K : PAdicLocalField p} {m : ℕ} {g : K.absGal} :
    g ∈ muFixer K m ↔ ∀ x : K.closure, x ^ m = 1 → g x = x := Iff.rfl

/-- **`muFixer` は正規**——`b⁻¹ x` も `m` 乗根だから。
★`K(ζ)` の固定部分群は一般には正規とは限らないので、共役作用を使う本経路では
`muFixer` の方を使わなければならない。 -/
instance muFixer_normal (K : PAdicLocalField p) (m : ℕ) : (muFixer K m).Normal := by
  constructor
  intro a ha b x hx
  have hb : (b⁻¹ x) ^ m = 1 := by rw [← map_pow, hx, map_one]
  show b (a (b⁻¹ x)) = x
  rw [ha _ hb]
  show (b * b⁻¹) x = x
  rw [mul_inv_cancel]; rfl

/-- 代数閉包には原始 `m` 乗根がある(標数 0)。 -/
theorem exists_isPrimitiveRoot_closure' (K : PAdicLocalField p) (m : ℕ) [NeZero m] :
    ∃ ζ : K.closure, IsPrimitiveRoot ζ m := by
  haveI : CharZero K.carrier :=
    charZero_of_injective_algebraMap (algebraMap ℚ_[p] K.carrier).injective
  haveI : NeZero ((m : ℕ) : K.carrier) := NeZero.charZero
  exact HasEnoughRootsOfUnity.exists_primitiveRoot K.closure m

/-- `ζ` が原始 `m` 乗根なら、`ζ` を止める元は `m` 乗根をすべて止める。 -/
theorem fixingSubgroup_le_muFixer (K : PAdicLocalField p) {m : ℕ} [NeZero m] {ζ : K.closure}
    (hζ : IsPrimitiveRoot ζ m) :
    (IntermediateField.adjoin K.carrier ({ζ} : Set K.closure)).fixingSubgroup ≤ muFixer K m := by
  intro g hg x hx
  have hgζ : g ζ = ζ := by
    rw [IntermediateField.mem_fixingSubgroup_iff] at hg
    exact hg ζ (IntermediateField.subset_adjoin K.carrier _ rfl)
  obtain ⟨i, _, rfl⟩ := hζ.eq_pow_of_pow_eq_one hx
  rw [map_pow, hgζ]

/-- **`muFixer K m` は開**——`K(ζ)` は `K` 上有限次なのでその固定部分群は開であり、
`muFixer` はそれを含む。 -/
theorem isOpen_muFixer (K : PAdicLocalField p) (m : ℕ) [NeZero m] :
    IsOpen ((muFixer K m : Subgroup K.absGal) : Set K.absGal) := by
  obtain ⟨ζ, hζ⟩ := exists_isPrimitiveRoot_closure' K m
  have hint : IsIntegral K.carrier ζ := IsAlgebraic.isIntegral (Algebra.IsAlgebraic.isAlgebraic _)
  haveI : FiniteDimensional K.carrier
      (IntermediateField.adjoin K.carrier ({ζ} : Set K.closure)) :=
    IntermediateField.adjoin.finiteDimensional hint
  exact Subgroup.isOpen_mono (fixingSubgroup_le_muFixer K hζ)
    (IntermediateField.fixingSubgroup_isOpen _)

/-! ## 2. 「`μ_m ⊆ L_H`」との言い換え -/

/-- **`H ≤ muFixer K m` は「`μ_m ⊆ L_H`」と同じこと**。 -/
theorem le_muFixer_iff (K : PAdicLocalField p) (m : ℕ) (S : Subgroup K.absGal) :
    S ≤ muFixer K m ↔ ∀ x : K.closure, x ^ m = 1 → x ∈ IntermediateField.fixedField S := by
  constructor
  · intro h x hx
    rw [IntermediateField.mem_fixedField_iff]
    exact fun f hf => h hf x hx
  · intro h g hg x hx
    exact (IntermediateField.mem_fixedField_iff _ _).mp (h x hx) g hg

/-- `H ≤ muFixer K m` なら固定体の中に原始 `m` 乗根がある。 -/
theorem exists_isPrimitiveRoot_fixedField_of_le_muFixer (K : PAdicLocalField p) {m : ℕ} [NeZero m]
    {S : Subgroup K.absGal} (hS : IsOpen (S : Set K.absGal)) (hle : S ≤ muFixer K m) :
    ∃ η : (fixedFieldLocalField K S hS).carrier, IsPrimitiveRoot η m := by
  obtain ⟨ζ, hζ⟩ := exists_isPrimitiveRoot_closure' K m
  have hmem : ζ ∈ IntermediateField.fixedField S :=
    (le_muFixer_iff K m S).mp hle ζ hζ.pow_eq_one
  exact ⟨⟨ζ, hmem⟩, IsPrimitiveRoot.of_map_of_injective
    (f := (IntermediateField.fixedField S).subtype) hζ (fun a b h => Subtype.ext h)⟩

/-! ## 3. ★Λ10 の本体——両側で `μ_m` を止める開正規部分群 -/

/-- **★★★★★★★★★★`μ_m ⊆ L_H` と `μ_m ⊆ L_{α(H)}` を同時に満たす
開正規部分群 `H` の構成**。

`H := muFixer K m ⊓ α⁻¹(muFixer K' m)`。開正規 ⊓ 開正規なので開正規。
★`H` は `∃` の内側に閉じ込められている——結論の statement には現れない。 -/
theorem exists_openNormal_le_muFixer_both (K K' : PAdicLocalField p)
    (α : ContinuousMulEquiv K.absGal K'.absGal) (m : ℕ) [NeZero m] :
    ∃ H : Subgroup K.absGal, H.Normal ∧ IsOpen (H : Set K.absGal) ∧
      H ≤ muFixer K m ∧ (H.map α.toMulEquiv.toMonoidHom) ≤ muFixer K' m := by
  haveI hBn : (Subgroup.comap α.toMulEquiv.toMonoidHom (muFixer K' m)).Normal :=
    Subgroup.Normal.comap (muFixer_normal K' m) _
  refine ⟨muFixer K m ⊓ Subgroup.comap α.toMulEquiv.toMonoidHom (muFixer K' m),
    Subgroup.normal_inf_normal _ _, ?_, inf_le_left, ?_⟩
  · rw [Subgroup.coe_inf]
    exact (isOpen_muFixer K m).inter
      (isOpen_comap_of_continuousMulEquiv α _ (isOpen_muFixer K' m))
  · rintro _ ⟨x, hx, rfl⟩
    exact hx.2

/-- **Λ10 を体の言葉で述べたもの**——`L_H` と `L_{α(H)}` が同時に原始 `m` 乗根を持つ
開正規部分群 `H` が存在する。 -/
theorem exists_openNormal_isPrimitiveRoot_both (K K' : PAdicLocalField p)
    (α : ContinuousMulEquiv K.absGal K'.absGal) (m : ℕ) [NeZero m] :
    ∃ (H : Subgroup K.absGal) (_ : H.Normal) (hH : IsOpen (H : Set K.absGal))
      (hH' : IsOpen ((H.map α.toMulEquiv.toMonoidHom : Subgroup K'.absGal) : Set K'.absGal)),
      (∃ η : (fixedFieldLocalField K H hH).carrier, IsPrimitiveRoot η m) ∧
      (∃ η' : (fixedFieldLocalField K' (H.map α.toMulEquiv.toMonoidHom) hH').carrier,
        IsPrimitiveRoot η' m) := by
  obtain ⟨H, hHn, hHopen, hHle, hH'le⟩ := exists_openNormal_le_muFixer_both K K' α m
  have hH'open : IsOpen ((H.map α.toMulEquiv.toMonoidHom : Subgroup K'.absGal) : Set K'.absGal) :=
    isOpen_map_of_continuousMulEquiv α H hHopen
  exact ⟨H, hHn, hHopen, hH'open,
    exists_isPrimitiveRoot_fixedField_of_le_muFixer K hHopen hHle,
    exists_isPrimitiveRoot_fixedField_of_le_muFixer K' hH'open hH'le⟩

/-! ## 4. 群論の仕上げ——同型類ではなく「作用の仕方」が指標を決める -/

/-- `ℤ/m` の元は、`Multiplicative (ℤ/m)` への冪作用が一致すれば等しい。 -/
theorem zmod_eq_of_multiplicative_pow_eq {m : ℕ} [NeZero m] {a b : ZMod m}
    (h : ∀ z : Multiplicative (ZMod m), z ^ a.val = z ^ b.val) : a = b := by
  have h1 := h (Multiplicative.ofAdd (1 : ZMod m))
  rw [← ofAdd_nsmul, ← ofAdd_nsmul] at h1
  have h2 : (a.val : ℕ) • (1 : ZMod m) = (b.val : ℕ) • (1 : ZMod m) :=
    Multiplicative.ofAdd.injective h1
  simp only [nsmul_eq_mul, mul_one] at h2
  rwa [ZMod.natCast_zmod_val, ZMod.natCast_zmod_val] at h2

/-- **★同型類の不定性は指標に伝わらない**。

`e₁ : Λ ≃* M`、`e₂ : Λ' ≃* M`、`T : Λ ≃* Λ'` があり、`σ` と `σ'` が `T` で対応し、
`e₁` が `σ` を `c` 乗に、`e₂` が `σ'` を `c'` 乗に写すなら、`M` 上で `c` 乗と `c'` 乗は一致する。

証明は `Φ := e₂ ∘ T ∘ e₁⁻¹` が群同型であること 1 点による:
`Φ(y^c) = (Φ y)^{c'}` かつ `Φ(y^c) = (Φ y)^c`。
★`e₁`・`e₂` の取り方(生成元の選び方)は結論に効かない——これが
「核だけでは `χ` は決まらないが、円分子の**作用**なら決まる」の中身である。 -/
theorem pow_eq_pow_of_equivariant {M Λ Λ' : Type*} [Group M] [Group Λ] [Group Λ']
    (e₁ : Λ ≃* M) (e₂ : Λ' ≃* M) (T : Λ ≃* Λ') (σ : Λ ≃* Λ) (σ' : Λ' ≃* Λ') (c c' : ℕ)
    (hT : ∀ x, T (σ x) = σ' (T x))
    (h1 : ∀ x, e₁ (σ x) = (e₁ x) ^ c)
    (h2 : ∀ y, e₂ (σ' y) = (e₂ y) ^ c') :
    ∀ w : M, w ^ c = w ^ c' := by
  intro w
  set Φ : M ≃* M := (e₁.symm.trans T).trans e₂ with hΦ
  obtain ⟨z, rfl⟩ : ∃ z, Φ z = w := ⟨Φ.symm w, by simp⟩
  have hx := h1 (e₁.symm z)
  rw [MulEquiv.apply_symm_apply] at hx
  have hσ : σ (e₁.symm z) = e₁.symm (z ^ c) := by rw [← hx, MulEquiv.symm_apply_apply]
  have key : e₂ (T (σ (e₁.symm z))) = (e₂ (T (e₁.symm z))) ^ c' := by rw [hT, h2]
  rw [hσ] at key
  have hL : Φ (z ^ c) = e₂ (T (e₁.symm (z ^ c))) := rfl
  have hR : Φ z = e₂ (T (e₁.symm z)) := rfl
  rw [← hL, ← hR, map_pow] at key
  exact key

/-! ## 5. ★残った壁——名前の付いた 1 本の補題 -/

/-- **★★★★★★★★★★★★★★★これだけが残った壁**——局所類体論が与えるはずの 1 本。

「`μ_{p^n} ⊆ L_S` なる開正規 `S ⊴ Γ_F` について、円分子 `Λ_n(S) = tors_{p^n}(S^{ab})` は
巡回群 `ℤ/p^n` と同型で、`Γ_F` の共役作用は円分指標 `χ_{F,n}` 倍として働く。」

古典的には `S^{ab} ≅ (L_S^×)^∧`(Lubin–Tate 相互律)から、
`L_S^× ≅ ℤ × μ(L_S) × ℤ_p^d` の `p^n` 捩れが `μ_{p^n}(L_S) = μ_{p^n}` であることによる。
相互律の Galois 同変性が「作用が `χ` 倍」を与える。

★**この主張は 1 つの体 `F` の中で閉じている**——2 つの体を比べていない。
したがって `cyclotomicRecovery_of_torsionCyclotome` に循環はない。

★退化の自己検査はモジュール docstring を参照(`IsOpen S` も
`S ≤ muFixer F (p^n)` も落とすと**偽になる**)。 -/
def TorsionCyclotomeIsCyclotomic (p : ℕ) [Fact p.Prime] : Prop :=
  ∀ (F : PAdicLocalField p) (n : ℕ) (S : Subgroup F.absGal) (hS : S.Normal),
    IsOpen (S : Set F.absGal) → S ≤ muFixer F (p ^ n) →
    ∃ e : ↥(cyclotome ↥S (p ^ n)) ≃* Multiplicative (ZMod (p ^ n)),
      ∀ (g : F.absGal) (x : ↥(cyclotome ↥S (p ^ n))),
        e (@cyclotomeConj F.absGal _ _ _ S hS (p ^ n) g x)
          = (e x) ^ ((PadicInt.toZModPow n)
              ((cyclotomicCharacter F.closure p g.toRingEquiv : ℤ_[p]))).val

/-- **`IsOpen S` は落とせない**——`S = ⊥` は正規で `S ≤ muFixer F m` を満たすが、
`Λ_n(⊥)` は自明群なので `p^n > 1` のとき `ℤ/p^n` と同型になりえない。 -/
theorem not_nonempty_cyclotomeEquiv_bot (Γ : Type*) [Group Γ] [TopologicalSpace Γ]
    [IsTopologicalGroup Γ] {m : ℕ} (hm : 1 < m) :
    ¬ Nonempty (↥(cyclotome ↥(⊥ : Subgroup Γ) m) ≃* Multiplicative (ZMod m)) := by
  rintro ⟨e⟩
  haveI : Subsingleton (TopologicalAbelianization ↥(⊥ : Subgroup Γ)) := by
    constructor
    intro x y
    induction x using QuotientGroup.induction_on with
    | H a => induction y using QuotientGroup.induction_on with
      | H b => rw [Subsingleton.elim a b]
  haveI : Subsingleton (Multiplicative (ZMod m)) := e.toEquiv.symm.subsingleton
  have h01 : (0 : ZMod m) = 1 :=
    Multiplicative.ofAdd.injective
      (Subsingleton.elim (Multiplicative.ofAdd (0 : ZMod m)) (Multiplicative.ofAdd (1 : ZMod m)))
  haveI : NeZero m := ⟨by omega⟩
  have hv := congrArg ZMod.val h01
  rw [ZMod.val_zero, ZMod.val_one'' (by omega)] at hv
  exact absurd hv (by omega)

/-! ## 6. ★到達点——Proposition 1.1 の還元 -/

/-- **★★★★★★★★★★★★★★★★★★★★★★★★Proposition 1.1 は
`TorsionCyclotomeIsCyclotomic p` に還元される**。

すなわち「`μ_{p^n}` を含む開正規部分群の円分子が `Γ`-加群として `μ_{p^n}` である」
という**局所類体論の 1 補題**さえあれば、原典 Proposition 1.1
(円分指標が `Γ_K` から群論的に回復できる)が従う。

証明の要点は 3 つで、いずれも本経路の設計そのものである:

1. `H := muFixer K (p^n) ⊓ α⁻¹(muFixer K' (p^n))` を取ると、`H` も `α(H)` も
   `μ_{p^n}` を止める開正規部分群になる(`exists_openNormal_le_muFixer_both` と同じ構成)。
   ★`H` は canonical でなくてよい——必要なのは `α` で対応することだけ。
2. `α` が誘導する円分子の同型は共役作用と同変(Λ2 の
   `cyclotomeEquiv_subgroupMapCME_conj`)。
3. 同型の取り方(生成元の選択)は指標に影響しない(`pow_eq_pow_of_equivariant`)。

★接続先は `Found/PGC/CyclotomicRecovery.lean` の
`cyclotomicCharacterObject_recoverable_iff`(sorry 無し)。
`cyclotomicCharacterObject_transport_of_moduleEquiv` が局所 Tate 双対性を
入力に取るのに対し、本定理は**局所類体論**を入力に取る——同じ結論への 2 本目の橋である。 -/
theorem cyclotomicRecovery_of_torsionCyclotome (hLCFT : TorsionCyclotomeIsCyclotomic p) :
    (cyclotomicCharacterObject (p := p)).RecoverableFromAbsGal := by
  refine cyclotomicCharacterObject_recoverable_iff.mpr ?_
  intro K K' α g
  apply Units.ext
  refine PadicInt.ext_of_toZModPow.mp ?_
  intro n
  haveI : NeZero (p ^ n) := ⟨pow_ne_zero _ (Fact.out : p.Prime).ne_zero⟩
  set A : Subgroup K.absGal := muFixer K (p ^ n) with hA
  set B : Subgroup K.absGal :=
    Subgroup.comap α.toMulEquiv.toMonoidHom (muFixer K' (p ^ n)) with hB
  haveI hBn : B.Normal := Subgroup.Normal.comap (muFixer_normal K' (p ^ n)) _
  set H : Subgroup K.absGal := A ⊓ B with hHdef
  haveI hHn : H.Normal := Subgroup.normal_inf_normal A B
  have hAopen : IsOpen (A : Set K.absGal) := isOpen_muFixer K (p ^ n)
  have hBopen : IsOpen (B : Set K.absGal) :=
    isOpen_comap_of_continuousMulEquiv α _ (isOpen_muFixer K' (p ^ n))
  have hHopen : IsOpen (H : Set K.absGal) := by
    rw [hHdef, Subgroup.coe_inf]; exact hAopen.inter hBopen
  set H' : Subgroup K'.absGal := H.map α.toMulEquiv.toMonoidHom with hH'def
  haveI hH'n : H'.Normal := hHn.map _ α.toMulEquiv.surjective
  have hH'open : IsOpen (H' : Set K'.absGal) := isOpen_map_of_continuousMulEquiv α H hHopen
  have hHle : H ≤ muFixer K (p ^ n) := inf_le_left
  have hH'le : H' ≤ muFixer K' (p ^ n) := by
    rintro _ ⟨x, hx, rfl⟩
    exact hx.2
  obtain ⟨e₁, he₁⟩ := hLCFT K n H hHn hHopen hHle
  obtain ⟨e₂, he₂⟩ := hLCFT K' n H' hH'n hH'open hH'le
  have hkey := pow_eq_pow_of_equivariant e₁ e₂ (cyclotomeEquiv (subgroupMapCME α H) (p ^ n))
    (cyclotomeConj H (p ^ n) g) (@cyclotomeConj K'.absGal _ _ _ H' hH'n (p ^ n) (α g))
    _ _ (cyclotomeEquiv_subgroupMapCME_conj α H hH'n (p ^ n) g) (he₁ g) (he₂ (α g))
  exact (zmod_eq_of_multiplicative_pow_eq hkey).symm

end ABC3.Found.PGC
