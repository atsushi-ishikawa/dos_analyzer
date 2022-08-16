from ase import Atoms, Atom
from ase.calculators.vasp import VaspDos
import sys
import numpy as np
import math
import matplotlib.pylab as plt

class VaspDosPlus:
    def __init__(self, doscar=None):
        self.doscar = doscar

        self._system = None
        self._numpeaks = None
        self.db = None
        self._relative_to_fermi = None

        # --- conditions ---
        self.margin = 0.0
        # Margin for the occupied band. Efermi + margin is considered as occupied.
        # margin = 0.0 is not good.
        # Better to increase this value when including edge.
        self.normalize_height = True    # True is maybe better
        self.do_hilbert = False
        self.do_cohp = False
        self.geometry_information = True

        #self.sigma = [20, 40, 70]  # smearing width ... 30-60
        #self.sigma = [20, 40, 80]  # smearing width ... 30-60
        self.sigma = [4, 4, 4]
        # --- end conditions ---

        self.vaspdos = VaspDos(doscar=doscar)
        self.efermi = get_efermi_from_doscar(self.doscar)

        # limit analysis on occupied part only
        energy = self.vaspdos.energy
        self.energy = np.array(list(filter(lambda x: x <= self.efermi + self.margin, energy)))

        # finding natom
        with open(doscar, "r") as f:
            line1 = f.readline()
            self.natom = int(line1.split()[0])

    @property
    def relative_to_fermi(self):
        return self._relative_to_fermi

    @relative_to_fermi.setter
    def relative_to_fermi(self, relative_to_fermi=False):
        self._relative_to_fermi = relative_to_fermi
        if relative_to_fermi:
            self.energy -= self.efermi

    @property
    def numpeaks(self):
        return self._numpeaks

    @numpeaks.setter
    def numpeaks(self, number_of_peaks=1):
        self._numpeaks = number_of_peaks

    @property
    def system(self):
        return self._system

    @system.setter
    def system(self, systemname=None):
        self._system = systemname

    def set_database(self, jsonfile=None):
        from ase.db import connect
        self.db = connect(jsonfile)

    def get_descriptors(self, adsorbate=False, with_centers=True):
        """
        Get descriptors. Descriptors are peak positions, widths, and heights of s, p, and d-bands.
        The number of peaks should be specified by VaspDosPlus.numpeaks.

        Args:
            adsorbate: whether adsorbate or not
            with_centers: whether to include s, p, d-band centers or not
        Returns:
            descriptors (dict)
        """
        check = False
        sort_key = "position"  # "height" or "position"

        descriptors = {}
        descriptors = self.add_system_to_dict(dict=descriptors)

        if self.do_cohp:
            cohpcar = self._system + "COHPCAR"

        print(" ----- system - {0:s} ---- ".format(self._system))

        if adsorbate:
            atom_range = range(0, self.natom)
            orbitals = {"s": 0, "p": 1}
        else:
            atom_range = range(48, 64)
            orbitals = {"s": 0, "p": 1, "d": 2}

        for orbital in orbitals.values():
            pdos = self.get_pdos(self.vaspdos, atom_range=atom_range, orbital=orbital)

            # smear the dos
            pdos = self.smear_dos(pdos, sigma=self.sigma[orbital])

            # find peaks to make guess for Gaussian fit
            peaks = self.findpeak(pdos)

            width = 1.0 * (1 / self.sigma[orbital]**0.5)  # guess
            params = []
            for idx in peaks:
                params.append(self.energy[idx])
                params.append(pdos[idx])
                params.append(width)

            # Try gaussian fit. If fails, return blanked list
            try:
                orb_name = self.from_012_to_spd(orbital)
                params, rss, r2 = gaussian_fit(np.array(self.energy), pdos, params)
                peaks = sort_peaks(params, key=sort_key)
                print("found {0:>2d} peaks -- {1:s} orbital R^2 = {2:>5.3f}".format(len(peaks), orb_name, r2))
            except:
                r2 = 0.0
                peaks = [(0, 0, 0) for _ in range(self._numpeaks)]

            # discard if R^2 is too low
            if r2 < 0.90:
                print("fitting failed: R^2 ({:>5.3f}) is too low ... quit".format(r2))
                peaks = [(0, 0, 0) for _ in range(self._numpeaks)]

            # sort peaks and limit to numpeaks
            peaks = peaks[0:self._numpeaks]
            tmp = []
            for peak in peaks:
                tmp.append(peak[0])
                tmp.append(peak[1])
                tmp.append(peak[2])

            peaks = sort_peaks(tmp, key=sort_key)

            if self._relative_to_fermi:
                occ_peaks = list(filter(lambda x: x[0] <  self.margin, peaks))
                vir_peaks = list(filter(lambda x: x[0] >= self.margin, peaks))
            else:
                occ_peaks = list(filter(lambda x: x[0] <  self.efermi + self.margin, peaks))
                vir_peaks = list(filter(lambda x: x[0] >= self.efermi + self.margin, peaks))

            # zero padding upto numpeaks
            occ_peaks = occ_peaks + [(0, 0, 0)] * (self._numpeaks - len(occ_peaks))
            vir_peaks = vir_peaks + [(0, 0, 0)] * (self._numpeaks - len(vir_peaks))

            if self.do_cohp:
                cohp_pos_peak, cohp_pos_center, cohp_neg_peak, cohp_neg_center = cohp_analysis(cohpcar)

            # occupied and virtual
            position_occ = []
            height_occ   = []
            width_occ    = []
            position_vir = []
            height_vir   = []
            width_vir    = []

            for peak in occ_peaks:
                position_occ.append(peak[0])
                height_occ.append(peak[1])
                width_occ.append(peak[2])

            for peak in vir_peaks:
                position_vir.append(peak[0])
                height_vir.append(peak[1])
                width_vir.append(peak[2])

            # if you want to check by eye
            if check:
                self.plot_dos(pdos, params=params, position_occ=position_occ)

            # write to database
            orb_name = self.from_012_to_spd(orbital)

            for ipeak in range(self._numpeaks):
                descriptors.update({orb_name + "_position_occ_" + str(ipeak): position_occ[ipeak]})
                descriptors.update({orb_name + "_height_occ_" + str(ipeak): height_occ[ipeak]})
                descriptors.update({orb_name + "_width_occ_" + str(ipeak): width_occ[ipeak]})

            # get s, p, d-center and higher moments
            if with_centers:
                center = self.get_moments(pdos, order=1)
                descriptors.update({orb_name + "_center": center})
            # higher moments ... not effective
            #second = self.get_moments(pdos, order=2)
            #descriptors.update({orb_name + "_second": second})
            #third = self.get_moments(pdos, order=3)
            #descriptors.update({orb_name + "_third": third})
            #forth = self.get_moments(pdos, order=3)
            #descriptors.update({orb_name + "_forth": forth})

        # end loop for orbitals

        if not self._relative_to_fermi:
            descriptors.update({"e_fermi": self.efermi})

        if self.do_cohp:
            descriptors.update({"cohp_pos_center": cohp_pos_center})
            descriptors.update({"cohp_neg_center": cohp_neg_center})
            descriptors.update({"cohp_pos_peak": cohp_pos_peak})
            descriptors.update({"cohp_neg_peak": cohp_neg_peak})

        # take upper edge of TOTAL DOS by inverse Hilbert transform
        if self.do_hilbert:
            tdos = self.vaspdos.dos
            upper_edge, lower_edge = self.get_dos_edge(tdos)
            if self._relative_to_fermi:
                upper_edge -= self.efermi
                lower_edge -= self.efermi

            descriptors.update({"upper_edge": upper_edge})
            descriptors.update({"lower_edge": lower_edge})

        return descriptors

    def add_system_to_dict(self, dict=None):
        """
        Add surface composition to the dict.

        Args:
            dict:
        Returns:
            dict:
        """
        tmp = {"system": self._system}
        dict.update(tmp)
        return dict

    def cohp_analysis(self, cohpcar):
        ene2, cohp = read_cohp(cohpcar)
        ene2 = np.array(ene2) + self.efermi

        cohp_pos = np.array(list(map(lambda x: x if x > 0.0 else 0.0, cohp)))
        cohp_neg = np.array(list(map(lambda x: x if x < 0.0 else 0.0, cohp)))

        cohp_pos = self.smear_dos(cohp_pos, sigma=self.sigma)
        cohp_neg = self.smear_dos(cohp_neg, sigma=self.sigma)

        # find peak for pos
        peaks = self.findpeak(cohp_pos)
        maxind = peaks[np.argmax(cohp_pos[peaks])]  # highest peak
        cohp_pos_peak = ene2[maxind]

        # find peak for neg
        peaks = self.findpeak(cohp_neg)
        minind = peaks[np.argmin(cohp_neg[peaks])]  # highest peak
        cohp_neg_peak = ene2[minind]

        cohp_pos_center = get_distribution_moment(cohp_pos, order=1)
        cohp_neg_center = get_distribution_moment(cohp_neg, order=1)

        return cohp_pos_peak, cohp_pos_center, cohp_neg_peak, cohp_neg_center

    def plot_dos(self, dos, params=None, position_occ=None):
        """
        Plot the density of state.

        Args:
            dos: density of state (numpy array)
            params:
            position_occ: list of occupied peaks
        """
        import seaborn as sb
        from scipy import fftpack

        fit = fit_func(self.energy, *params)

        if self.do_hilbert:
            ih = fftpack.ihilbert(pdos)

        sb.set(context='notebook', style='darkgrid', palette='deep',font='sans-serif', font_scale=1,
               color_codes=False, rc=None)
        plt.plot(self.energy, fit, label="fitted")
        plt.plot(self.energy, dos, label="original")
        plt.vlines(x=position_occ, ymin=0, ymax=np.max(dos), linestyle="dashed", linewidth=0.5)

        if self.do_hilbert:
            plt.plot(self.energy, ih, label="inverse Hilbert")

        if self._relative_to_fermi:
            plt.xlim([min(self.energy), 0])
        else:
            plt.xlim([min(self.energy), max(self.energy) + self.margin])

        plt.legend()
        plt.show()

        if self.do_cohp:
            cohp_smeared = self.smear_dos(cohp, sigma=self.sigma)
            plt.plot(ene2, cohp_smeared, label="COHP")
            plt.legend()
            plt.show()

        return None

    def get_moments(self, pdos=None, order=None):
        """
        Get moments.

        Args:
            pdos:
            order:
        Returns:
            order-th moment value (float)
        """
        from ase.dft import get_distribution_moment
        return get_distribution_moment(self.energy, pdos, order=order)

    def get_adsorption_site(self, atoms=None, adsorbing_element="C"):
        """
        Get coordinating atom.

        Args:
            atoms:
            adsorbing_element:
        Returns:
            coord_num (int): the atom index of the adsorption site.
        """
        coord_ind = []
        adsorbate_ind = atoms.get_chemical_symbols().index(adsorbing_element)
        li = [atoms.get_distance(adsorbate_ind, i) for i in range(len(atoms))]
        li = li[:adsorbate_ind]  # limit on surface
        li = list(map(lambda x: 1.0e3 if x < 1.0e-3 else x, li))

        coord_ind.append(np.argmin(li))
        coord_num = atoms[coord_ind].get_atomic_numbers()[0]  # atomic number of coordinating atom

        return coord_num

    def get_pdos(self, vaspdos, atom_range=None, orbital=None):
        """
        Get projected DOS.

        Args:
            vaspdos: VaspDos object
            atom_range:
            orbital:
        Returns:
            pdos: projected DOS.
        """
        pdos = np.zeros(len(self.energy))
        for i in atom_range:
            pdos += vaspdos.site_dos(i, orbital)[:len(self.energy)]

        pdos = pdos[:len(self.energy)]

        if self.normalize_height:
            pdos = pdos / np.max(pdos)

        return pdos

    def gaussian(self, x, x0, a, b):
        """
        y = a*exp(-b*(x-x0)**2)

        Args:
            x:
            x0:
            a:
            b:
        Returns:
            y
        """
        x = np.array(x)
        y = np.exp(-b * (x - x0)**2)
        y = a * y
        return y

    def smear_dos(self, dos, sigma=5.0):
        """
        Get smeared DOS.

        Args:
            dos:
            sigma:
        Returns:
            dos:
        """
        x = self.energy  # note: occpied part only
        y = dos
        smeared = np.zeros(len(x))

        for i, j in enumerate(x):
            x0 = x[i]
            a = y[i]
            smeared += self.gaussian(x, x0, a=a, b=sigma)

        return smeared

    def findpeak(self, y):
        """
        Find peak.

        Args:
            y:
        Returns:
            index
        """
        import peakutils
        indexes = peakutils.indexes(y, thres=0.1, min_dist=1)
        return indexes

    def from_012_to_spd(self, val):
        """
        Convert 012 to spd.

        Args:
            val:
        Returns:
            s or p or d (string)
        """
        spd = {"s": 0, "p": 1, "d": 2}
        keys = [k for k, v in spd.items() if v == val]
        if keys:
            return keys[0]

    def get_dos_edge(self, tdos):
        """
        Get the edge of the DOS, detected by the peak in its inverse Hilbert transform.
        Note: assuming Total DOS, but PDOS might be OK.

        Args:
            tdos:
        Returns:
            upper_edge, lower_edge (float):
        """
        from scipy import fftpack

        tdos = tdos[:len(self.energy)]
        tdos = self.smear_dos(tdos)

        # take upper edge by inverse Hilbert transform
        ih = fftpack.ihilbert(tdos)
        ihpeaks = self.findpeak(abs(ih))

        upper_edge = np.zeros(self._numpeaks)
        lower_edge = np.zeros(self._numpeaks)

        for peak in ihpeaks:
            if ih[peak] > 0.0 and ih[peak] > 0.8 * max(ih):  # just choose large peak, positive part
                upper_edge = np.insert(upper_edge, 0, self.energy[peak])
            elif ih[peak] <= 0.0 and ih[peak] < 0.8 * min(ih):
                lower_edge = np.insert(lower_edge, 0, self.energy[peak])

        #upper_edge = upper_edge[0:self._numpeaks]
        #lower_edge = lower_edge[0:self._numpeaks]
        #upper_edge = upper_edge[::-1]
        #lower_edge = lower_edge[::-1]
        upper_edge = upper_edge[0]
        lower_edge = lower_edge[0]

        return upper_edge, lower_edge


def get_efermi_from_doscar(DOSCAR):
    """
    Get Fermi energy from DOSCAR file.

    Args:
        DOSCAR:
    Returns:
        efermi: Fermi energy
    """
    import linecache
    line = linecache.getline(DOSCAR, 6)
    line = line.strip()
    line = line.split()
    efermi = float(line[3])
    return efermi

def gaussian_fit(x, y, guess):
    """

    Args:
        x:
        y:
        guess:

    Returns:

    """
    from scipy.optimize import curve_fit

    # method: trf ... Trust region reflective method. Generally robust (default).
    #         lm  ... Levenberg-Marquardt. Most efficient for small sparse problem.
    # ftol, xtol, gtol: default is 1.0e-8. Used 1.0e-6 to reduced the drop-off DOSs.
    #
    tol = 1.0e-8   # OK
    popt, pcov = curve_fit(fit_func, x, y, p0=guess, method="trf", ftol=tol, xtol=tol, gtol=tol)
    #popt, pcov = curve_fit(fit_func, x, y, p0=guess, method="lm", ftol=tol, xtol=tol, gtol=tol)

    fit = fit_func(x, *popt)
    residual = y - fit
    rss = np.sum(residual**2)  # residual sum of squares
    tss = np.sum((y - np.mean(y))**2)  # total sum of squares
    r2 = 1 - (rss / tss)

    return popt, rss, r2

def fit_func(x, *params):
    """

    Args:
        x:
        *params:
    Returns:

    """
    y = np.zeros_like(x)
    for i in range(0, len(params), 3):
        ctr = params[i]  # position
        amp = params[i + 1]  # height
        wid = params[i + 2]  # width (smaller is sharper)
        y = y + amp * np.exp(-((x - ctr) / wid) ** 2)
    return y

def sort_peaks(peaks, key="height"):
    """
    Sort peaks.
    Assuming peaks are stored in [position, height, width,  position, height, width,...]

    Args:
        peaks:
        key: height or position
    Returns:

    """
    dtype = [("position", float), ("height", float), ("width", float)]

    newpeaks = np.array([], dtype=dtype)

    for i in range(0, len(peaks), 3):
        peak = np.array((peaks[i], peaks[i + 1], peaks[i + 2]), dtype=dtype)
        newpeaks = np.append(newpeaks, peak)

    newpeaks = np.sort(newpeaks, order=key)
    newpeaks = newpeaks[::-1]  # sort in descending order
    return newpeaks
